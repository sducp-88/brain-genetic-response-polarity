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
    variant = NULL,
    min_cells = 20L,
    workers = 2L
  )
  i <- 1L
  while (i <= length(x)) {
    key <- gsub("-", "_", sub("^--", "", x[[i]]))
    if (i == length(x)) stop("Missing value for --", key)
    out[[key]] <- x[[i + 1L]]
    i <- i + 2L
  }
  for (key in c("input_dir", "output_dir", "variant")) {
    if (is.null(out[[key]])) stop("--", gsub("_", "-", key), " is required")
  }
  allowed <- c("exclude_age89plus", "omit_nonad", "omit_log_ncells")
  if (!out$variant %in% allowed) {
    stop("--variant must be one of: ", paste(allowed, collapse = ", "))
  }
  out$min_cells <- as.integer(out$min_cells)
  out$workers <- as.integer(out$workers)
  out
}

zscore <- function(x) {
  value <- as.numeric(scale(as.numeric(x)))
  if (any(!is.finite(value))) stop("Cannot standardize a constant covariate")
  value
}

age_numeric <- function(x) {
  as.numeric(sub("\\+$", "", as.character(x)))
}

ancestry_group <- function(x) {
  value <- as.character(x)
  factor(
    ifelse(value %in% c("AFR", "AMR", "EUR"), value, "Other"),
    levels = c("AFR", "AMR", "EUR", "Other")
  )
}

non_ad_diagnosis <- function(x) {
  value <- trimws(as.character(x))
  value[is.na(value)] <- ""
  as.integer(value != "" & value != "AD")
}

write_csv_atomic <- function(x, path) {
  temporary <- paste0(path, ".tmp")
  write.csv(x, temporary, row.names = FALSE, quote = TRUE, na = "")
  if (file.exists(path)) file.remove(path)
  if (!file.rename(temporary, path)) stop("Failed atomic rename: ", path)
}

prepare_metadata <- function(cd, phenotype_name, min_cells, variant) {
  frame <- as.data.frame(cd)
  frame$sample_id <- rownames(frame)
  frame$n_cells_numeric <- as.numeric(frame$n_cells)
  frame$age_raw <- trimws(as.character(frame$Age))
  frame$age_numeric <- age_numeric(frame$Age)
  frame$pmi_numeric <- as.numeric(frame$PMI)
  frame$sex_model <- factor(as.character(frame$Sex), levels = c("F", "M"))
  frame$ancestry_model <- ancestry_group(frame$Ancestry)
  frame$nonAD_dx <- non_ad_diagnosis(frame$Diagnosis)
  frame$phenotype_raw <- if (phenotype_name == "CERAD") {
    as.numeric(frame$CERAD)
  } else if (phenotype_name == "Braak") {
    as.numeric(frame$Braak)
  } else {
    stop("Unsupported phenotype")
  }
  required <- c(
    "n_cells_numeric", "age_numeric", "pmi_numeric", "sex_model",
    "ancestry_model", "phenotype_raw"
  )
  keep <- complete.cases(frame[, required, drop = FALSE]) &
    frame$n_cells_numeric >= min_cells
  if (variant == "exclude_age89plus") {
    keep <- keep & frame$age_raw != "89+"
  }
  frame <- droplevels(frame[keep, , drop = FALSE])
  if (nrow(frame) < 40L) stop("Too few complete eligible donors")
  frame$phenotype_value <- zscore(frame$phenotype_raw)
  frame$age_z <- zscore(frame$age_numeric)
  frame$pmi_z <- zscore(frame$pmi_numeric)
  frame$log10_n_cells_z <- zscore(log10(frame$n_cells_numeric))
  frame
}

build_design <- function(metadata, variant) {
  terms <- c(
    "phenotype_value", "age_z", "sex_model", "pmi_z", "ancestry_model"
  )
  if (variant != "omit_log_ncells") {
    terms <- c(terms, "log10_n_cells_z")
  }
  if (variant != "omit_nonad") {
    terms <- c(terms, "nonAD_dx")
  }
  formula <- as.formula(paste("~", paste(terms, collapse = " + ")))
  design <- model.matrix(formula, data = metadata)
  if (qr(design)$rank != ncol(design)) stop("Design is not full rank")
  design
}

fit_one <- function(
    counts, gene_metadata, metadata, phenotype_name, class_name, variant,
    output_dir
) {
  design <- build_design(metadata, variant)
  sample_columns <- match(metadata$sample_id, colnames(counts))
  if (anyNA(sample_columns)) stop("Sample identifiers do not match")
  model_counts <- counts[, sample_columns, drop = FALSE]
  storage.mode(model_counts@x) <- "double"

  dge <- DGEList(counts = model_counts)
  minimum_samples <- max(10L, ceiling(0.10 * ncol(dge)))
  keep_gene <- rowSums(cpm(dge) >= 1) >= minimum_samples
  if (sum(keep_gene) < 1000L) stop("Implausibly few retained genes")
  dge <- dge[keep_gene, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  voom_object <- voom(
    dge, design = design, plot = FALSE, normalize.method = "none"
  )
  fit <- eBayes(lmFit(voom_object, design), robust = TRUE)
  coefficient <- which(colnames(design) == "phenotype_value")
  table <- topTable(
    fit, coef = coefficient, number = Inf, sort.by = "none",
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
      min_cells = 20L,
      variant = variant,
      stringsAsFactors = FALSE
    )
  )
  result_path <- file.path(
    output_dir, sprintf("%s__%s__gene_results.csv", class_name, phenotype_name)
  )
  write_csv_atomic(result, result_path)
  data.frame(
    class = class_name,
    phenotype = phenotype_name,
    status = "complete",
    variant = variant,
    voom_method = "voom",
    n_donors = nrow(metadata),
    age_min = min(metadata$age_numeric),
    age_max = max(metadata$age_numeric),
    genes_tested = nrow(result),
    design_rank = qr(design)$rank,
    design_columns = ncol(design),
    design_condition_number = kappa(design, exact = TRUE),
    fdr05_genes = sum(result$FDR_within_class < 0.05, na.rm = TRUE),
    result_file = basename(result_path),
    stringsAsFactors = FALSE
  )
}

settings <- parse_args(commandArgs(trailingOnly = TRUE))
input_dir <- normalizePath(settings$input_dir, mustWork = TRUE)
output_dir <- normalizePath(settings$output_dir, mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

classes <- c("Astro", "EN", "Immune", "IN", "Oligo", "OPC")
phenotypes <- c("CERAD", "Braak")
manifest_rows <- list()

for (class_name in classes) {
  input_file <- file.path(
    input_dir, sprintf("RADC_%s_pseudobulk_counts.h5ad", class_name)
  )
  message("Reading ", input_file)
  sce <- readH5AD(input_file, use_hdf5 = FALSE, reader = "R")
  counts <- assay(sce, "X")
  if (!inherits(counts, "sparseMatrix")) counts <- as(counts, "dgCMatrix")
  rd <- as.data.frame(rowData(sce))
  gene_metadata <- data.frame(
    gene_id = rownames(sce),
    gene_name = if ("gene_name" %in% colnames(rd)) {
      as.character(rd$gene_name)
    } else {
      rownames(sce)
    },
    stringsAsFactors = FALSE
  )
  rows <- parallel::mclapply(
    phenotypes,
    function(phenotype_name) {
      message(
        "Fitting ", settings$variant, " / ", class_name, " / ", phenotype_name
      )
      metadata <- prepare_metadata(
        colData(sce), phenotype_name, settings$min_cells, settings$variant
      )
      fit_one(
        counts, gene_metadata, metadata, phenotype_name, class_name,
        settings$variant, output_dir
      )
    },
    mc.cores = min(settings$workers, length(phenotypes)),
    mc.preschedule = FALSE
  )
  failed <- vapply(rows, inherits, logical(1), what = "try-error")
  if (any(failed)) stop(paste(as.character(rows[failed]), collapse = " | "))
  manifest_rows <- c(manifest_rows, rows)
  rm(sce, counts)
  gc()
}

manifest <- do.call(rbind, manifest_rows)
write_csv_atomic(manifest, file.path(output_dir, "model_manifest.csv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("COMPLETE ", settings$variant, ": ", output_dir)
