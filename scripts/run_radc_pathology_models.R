#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(zellkonverter)
  library(edgeR)
  library(limma)
  library(Matrix)
})

parse_args <- function(x) {
  out <- list(
    input_dir = NULL,
    output_dir = NULL,
    min_cells = 20L,
    cpm_threshold = 1,
    min_gene_samples = 10L,
    min_gene_fraction = 0.10,
    workers = 2L,
    classes = "Astro,EN,Endo,IN,Immune,Mural,OPC,Oligo",
    phenotypes = "CERAD,Braak,Dementia,AD_status",
    quality_weights = FALSE,
    overwrite = FALSE
  )
  i <- 1L
  while (i <= length(x)) {
    key <- sub("^--", "", x[[i]])
    if (key %in% c("overwrite", "quality-weights")) {
      out[[gsub("-", "_", key)]] <- TRUE
      i <- i + 1L
    } else {
      if (i == length(x)) stop("Missing value for --", key)
      value <- x[[i + 1L]]
      key <- gsub("-", "_", key)
      out[[key]] <- value
      i <- i + 2L
    }
  }
  if (is.null(out$input_dir) || is.null(out$output_dir)) {
    stop("--input-dir and --output-dir are required")
  }
  out$min_cells <- as.integer(out$min_cells)
  out$cpm_threshold <- as.numeric(out$cpm_threshold)
  out$min_gene_samples <- as.integer(out$min_gene_samples)
  out$min_gene_fraction <- as.numeric(out$min_gene_fraction)
  out$workers <- as.integer(out$workers)
  out$classes <- trimws(strsplit(out$classes, ",", fixed = TRUE)[[1]])
  allowed_classes <- c(
    "Astro", "EN", "Endo", "IN", "Immune", "Mural", "OPC", "Oligo"
  )
  if (
    !length(out$classes) ||
    any(!nzchar(out$classes)) ||
    any(!out$classes %in% allowed_classes) ||
    anyDuplicated(out$classes)
  ) {
    stop(
      "--classes must be a unique comma-separated subset of: ",
      paste(allowed_classes, collapse = ",")
    )
  }
  out$phenotypes <- trimws(strsplit(out$phenotypes, ",", fixed = TRUE)[[1]])
  allowed_phenotypes <- c("CERAD", "Braak", "Dementia", "AD_status")
  if (
    !length(out$phenotypes) ||
    any(!nzchar(out$phenotypes)) ||
    any(!out$phenotypes %in% allowed_phenotypes) ||
    anyDuplicated(out$phenotypes)
  ) {
    stop(
      "--phenotypes must be a unique comma-separated subset of: ",
      paste(allowed_phenotypes, collapse = ",")
    )
  }
  out
}

zscore <- function(x) {
  x <- as.numeric(x)
  value <- as.numeric(scale(x))
  if (any(!is.finite(value))) stop("Cannot standardize a constant covariate")
  value
}

age_numeric <- function(x) {
  as.numeric(sub("\\+$", "", as.character(x)))
}

ancestry_group <- function(x) {
  x <- as.character(x)
  out <- ifelse(x %in% c("AFR", "AMR", "EUR"), x, "Other")
  factor(out, levels = c("AFR", "AMR", "EUR", "Other"))
}

non_ad_diagnosis <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  as.integer(x != "" & x != "AD")
}

write_csv_atomic <- function(x, path) {
  temporary <- paste0(path, ".tmp")
  write.csv(x, temporary, row.names = FALSE, quote = TRUE, na = "")
  if (file.exists(path)) file.remove(path)
  if (!file.rename(temporary, path)) stop("Failed atomic rename: ", path)
}

extract_gene_metadata <- function(sce) {
  rd <- as.data.frame(rowData(sce))
  gene_id <- rownames(sce)
  gene_name <- if ("gene_name" %in% colnames(rd)) {
    as.character(rd$gene_name)
  } else if ("feature_name" %in% colnames(rd)) {
    as.character(rd$feature_name)
  } else {
    gene_id
  }
  data.frame(gene_id = gene_id, gene_name = gene_name, stringsAsFactors = FALSE)
}

prepare_metadata <- function(cd, phenotype_name, min_cells) {
  frame <- as.data.frame(cd)
  frame$sample_id <- rownames(frame)
  frame$n_cells_numeric <- as.numeric(frame$n_cells)
  frame$age_numeric <- age_numeric(frame$Age)
  frame$pmi_numeric <- as.numeric(frame$PMI)
  frame$sex_model <- factor(as.character(frame$Sex), levels = c("F", "M"))
  frame$ancestry_model <- ancestry_group(frame$Ancestry)
  frame$nonAD_dx <- non_ad_diagnosis(frame$Diagnosis)
  frame$phenotype_raw <- switch(
    phenotype_name,
    CERAD = as.numeric(frame$CERAD),
    Braak = as.numeric(frame$Braak),
    Dementia = as.numeric(frame$Dementia),
    AD_status = ifelse(
      as.character(frame$AD_status) == "Yes", 1,
      ifelse(as.character(frame$AD_status) == "No", 0, NA_real_)
    ),
    stop("Unknown phenotype: ", phenotype_name)
  )
  required <- c(
    "n_cells_numeric", "age_numeric", "pmi_numeric", "sex_model",
    "ancestry_model", "phenotype_raw"
  )
  complete <- complete.cases(frame[, required, drop = FALSE])
  keep <- complete & frame$n_cells_numeric >= min_cells
  frame <- droplevels(frame[keep, , drop = FALSE])
  if (nrow(frame) < 40L) stop("Too few complete eligible samples")
  if (length(unique(frame$phenotype_raw)) < 2L) stop("Constant phenotype")
  frame$phenotype_value <- if (phenotype_name %in% c("CERAD", "Braak")) {
    zscore(frame$phenotype_raw)
  } else {
    frame$phenotype_raw
  }
  frame$age_z <- zscore(frame$age_numeric)
  frame$pmi_z <- zscore(frame$pmi_numeric)
  frame$log10_n_cells_z <- zscore(log10(frame$n_cells_numeric))
  frame
}

fit_one_model <- function(
    counts, gene_metadata, metadata, phenotype_name, class_name, settings,
    output_dir
) {
  model_formula <- ~ phenotype_value + age_z + sex_model + pmi_z +
    ancestry_model + log10_n_cells_z + nonAD_dx
  design <- model.matrix(model_formula, data = metadata)
  if (qr(design)$rank != ncol(design)) {
    stop("Design is not full rank for ", class_name, " / ", phenotype_name)
  }
  condition_number <- kappa(design, exact = TRUE)
  sample_columns <- match(metadata$sample_id, colnames(counts))
  if (anyNA(sample_columns)) stop("Sample identifiers do not match count columns")
  model_counts <- counts[, sample_columns, drop = FALSE]
  storage.mode(model_counts@x) <- "double"

  dge <- DGEList(counts = model_counts)
  minimum_samples <- max(
    settings$min_gene_samples,
    ceiling(settings$min_gene_fraction * ncol(dge))
  )
  keep_gene <- rowSums(cpm(dge) >= settings$cpm_threshold) >= minimum_samples
  if (sum(keep_gene) < 1000L) stop("Implausibly few retained genes")
  dge <- dge[keep_gene, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")

  voom_object <- if (settings$quality_weights) {
    voomWithQualityWeights(
      dge,
      design = design,
      plot = FALSE,
      normalize.method = "none"
    )
  } else {
    voom(
      dge,
      design = design,
      plot = FALSE,
      normalize.method = "none"
    )
  }
  fit <- lmFit(voom_object, design)
  fit <- eBayes(fit, robust = TRUE)
  coefficient <- which(colnames(design) == "phenotype_value")
  if (length(coefficient) != 1L) stop("Phenotype coefficient not found")
  table <- topTable(
    fit,
    coef = coefficient,
    number = Inf,
    sort.by = "none",
    adjust.method = "BH"
  )
  moderated_t <- as.numeric(fit$t[, coefficient])
  moderated_se <- ifelse(
    is.finite(moderated_t) & moderated_t != 0,
    abs(as.numeric(fit$coefficients[, coefficient]) / moderated_t),
    NA_real_
  )
  result <- cbind(
    gene_metadata[keep_gene, , drop = FALSE],
    data.frame(
      class = class_name,
      phenotype = phenotype_name,
      beta_D = as.numeric(fit$coefficients[, coefficient]),
      SE_D_moderated = moderated_se,
      t = moderated_t,
      AveExpr = table$AveExpr,
      P.Value = table$P.Value,
      FDR_within_class = table$adj.P.Val,
      B = table$B,
      n_donors = nrow(metadata),
      min_cells = settings$min_cells,
      stringsAsFactors = FALSE
    )
  )

  sample_weights <- if (
    settings$quality_weights &&
    "sample.weights" %in% colnames(voom_object$targets)
  ) {
    as.numeric(voom_object$targets$sample.weights)
  } else {
    rep(1, nrow(metadata))
  }
  weights <- data.frame(
    sample_id = metadata$sample_id,
    donor_id = as.character(metadata$donor_id),
    class = class_name,
    phenotype = phenotype_name,
    phenotype_raw = metadata$phenotype_raw,
    n_cells = metadata$n_cells_numeric,
    sample_quality_weight = sample_weights,
    stringsAsFactors = FALSE
  )

  result_path <- file.path(
    output_dir, sprintf("%s__%s__gene_results.csv", class_name, phenotype_name)
  )
  weight_path <- file.path(
    output_dir, sprintf("%s__%s__sample_weights.csv", class_name, phenotype_name)
  )
  write_csv_atomic(result, result_path)
  write_csv_atomic(weights, weight_path)

  data.frame(
    class = class_name,
    phenotype = phenotype_name,
    status = "complete",
    voom_method = if (settings$quality_weights) {
      "voomWithQualityWeights"
    } else {
      "voom"
    },
    n_donors = nrow(metadata),
    phenotype_min = min(metadata$phenotype_raw),
    phenotype_max = max(metadata$phenotype_raw),
    phenotype_cases = if (phenotype_name %in% c("Dementia", "AD_status")) {
      sum(metadata$phenotype_raw == 1)
    } else {
      NA_integer_
    },
    phenotype_controls = if (phenotype_name %in% c("Dementia", "AD_status")) {
      sum(metadata$phenotype_raw == 0)
    } else {
      NA_integer_
    },
    genes_tested = nrow(result),
    design_rank = qr(design)$rank,
    design_columns = ncol(design),
    design_condition_number = condition_number,
    minimum_gene_samples = minimum_samples,
    sample_weight_min = if (all(is.na(sample_weights))) NA else min(sample_weights, na.rm = TRUE),
    sample_weight_median = if (all(is.na(sample_weights))) NA else median(sample_weights, na.rm = TRUE),
    sample_weight_max = if (all(is.na(sample_weights))) NA else max(sample_weights, na.rm = TRUE),
    fdr05_genes = sum(result$FDR_within_class < 0.05, na.rm = TRUE),
    fdr01_genes = sum(result$FDR_within_class < 0.01, na.rm = TRUE),
    positive_fdr05 = sum(
      result$FDR_within_class < 0.05 & result$beta_D > 0, na.rm = TRUE
    ),
    negative_fdr05 = sum(
      result$FDR_within_class < 0.05 & result$beta_D < 0, na.rm = TRUE
    ),
    result_file = basename(result_path),
    weight_file = basename(weight_path),
    stringsAsFactors = FALSE
  )
}

settings <- parse_args(commandArgs(trailingOnly = TRUE))
input_dir <- normalizePath(settings$input_dir, mustWork = TRUE)
output_dir <- normalizePath(settings$output_dir, mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
classes <- settings$classes
phenotypes <- settings$phenotypes

manifest_rows <- list()
all_result_paths <- unlist(lapply(classes, function(class_name) {
  file.path(
    output_dir,
    sprintf("%s__%s__gene_results.csv", class_name, phenotypes)
  )
}), use.names = FALSE)
for (class_name in classes) {
  input_file <- file.path(
    input_dir, sprintf("RADC_%s_pseudobulk_counts.h5ad", class_name)
  )
  message("Reading ", input_file)
  sce <- readH5AD(input_file, use_hdf5 = FALSE, reader = "R")
  if (!"X" %in% assayNames(sce)) stop("X assay missing: ", class_name)
  counts <- assay(sce, "X")
  if (!inherits(counts, "sparseMatrix")) counts <- as(counts, "dgCMatrix")
  if (nrow(counts) != 34176L) stop("Unexpected gene count: ", class_name)
  if (any(counts@x < 0) || any(abs(counts@x - round(counts@x)) > 1e-6)) {
    stop("Counts are not nonnegative integers: ", class_name)
  }
  gene_metadata <- extract_gene_metadata(sce)
  cd <- colData(sce)
  phenotypes_to_fit <- phenotypes[!vapply(phenotypes, function(phenotype_name) {
    result_file <- file.path(
      output_dir, sprintf("%s__%s__gene_results.csv", class_name, phenotype_name)
    )
    if (file.exists(result_file) && !settings$overwrite) {
      message("Skipping existing result ", basename(result_file))
      TRUE
    } else {
      FALSE
    }
  }, logical(1))]
  if (length(phenotypes_to_fit)) {
    class_rows <- parallel::mclapply(
      phenotypes_to_fit,
      function(phenotype_name) {
        message("Fitting ", class_name, " / ", phenotype_name)
        metadata <- prepare_metadata(cd, phenotype_name, settings$min_cells)
        fit_one_model(
          counts, gene_metadata, metadata, phenotype_name, class_name, settings,
          output_dir
        )
      },
      mc.cores = min(settings$workers, length(phenotypes_to_fit)),
      mc.preschedule = FALSE
    )
    failed <- vapply(class_rows, inherits, logical(1), what = "try-error")
    if (any(failed)) {
      stop(
        "Parallel model failure for ", class_name, ": ",
        paste(as.character(class_rows[failed]), collapse = " | ")
      )
    }
    manifest_rows <- c(manifest_rows, class_rows)
  }
  rm(sce, counts)
  gc()
}

if (length(manifest_rows)) {
  new_manifest <- do.call(rbind, manifest_rows)
  manifest_path <- file.path(output_dir, "model_manifest.csv")
  if (file.exists(manifest_path) && !settings$overwrite) {
    old_manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
    new_manifest <- rbind(old_manifest, new_manifest)
    new_manifest <- new_manifest[
      !duplicated(new_manifest[, c("class", "phenotype")], fromLast = TRUE), ,
      drop = FALSE
    ]
  }
  write_csv_atomic(new_manifest, manifest_path)
}

existing_results <- all_result_paths[file.exists(all_result_paths)]
if (length(existing_results) == length(classes) * length(phenotypes)) {
  combined <- do.call(rbind, lapply(existing_results, function(path) {
    read.csv(path, stringsAsFactors = FALSE)
  }))
  combined$FDR_global_phenotype <- ave(
    combined$P.Value,
    combined$phenotype,
    FUN = function(x) p.adjust(x, method = "BH")
  )
  write_csv_atomic(combined, file.path(output_dir, "all_gene_results.csv"))
}

sink(file.path(output_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()
message("RADC pathology model pipeline complete: ", output_dir)
