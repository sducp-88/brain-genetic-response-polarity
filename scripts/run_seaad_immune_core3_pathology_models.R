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
    input = NULL,
    output_dir = NULL,
    cpm_threshold = 1,
    min_gene_samples = 10L,
    min_gene_fraction = 0.10,
    cell_class = "Immune",
    phenotypes = "CERAD,Braak",
    quality_weights = FALSE,
    include_race = FALSE,
    overwrite = FALSE
  )
  i <- 1L
  while (i <= length(x)) {
    key <- sub("^--", "", x[[i]])
    if (key %in% c("quality-weights", "include-race", "overwrite")) {
      out[[gsub("-", "_", key)]] <- TRUE
      i <- i + 1L
    } else {
      if (i == length(x)) stop("Missing value for --", key)
      out[[gsub("-", "_", key)]] <- x[[i + 1L]]
      i <- i + 2L
    }
  }
  if (is.null(out$input) || is.null(out$output_dir)) {
    stop("--input and --output-dir are required")
  }
  out$cpm_threshold <- as.numeric(out$cpm_threshold)
  out$min_gene_samples <- as.integer(out$min_gene_samples)
  out$min_gene_fraction <- as.numeric(out$min_gene_fraction)
  allowed_classes <- c(
    "Astro", "EN", "Endo", "IN", "Immune", "Mural", "OPC", "Oligo"
  )
  if (!out$cell_class %in% allowed_classes) {
    stop(
      "--cell-class must be one of: ",
      paste(allowed_classes, collapse = ",")
    )
  }
  out$phenotypes <- trimws(strsplit(out$phenotypes, ",", fixed = TRUE)[[1]])
  if (
    !length(out$phenotypes) ||
    any(!out$phenotypes %in% c("CERAD", "Braak")) ||
    anyDuplicated(out$phenotypes)
  ) {
    stop("--phenotypes must be a unique subset of CERAD,Braak")
  }
  out
}

write_csv_atomic <- function(x, path) {
  temporary <- paste0(path, ".tmp")
  write.csv(x, temporary, row.names = FALSE, quote = TRUE, na = "")
  if (file.exists(path)) file.remove(path)
  if (!file.rename(temporary, path)) stop("Failed atomic rename: ", path)
}

zscore <- function(x) {
  result <- as.numeric(scale(as.numeric(x)))
  if (any(!is.finite(result))) stop("Cannot standardize covariate")
  result
}

prepare_metadata <- function(cd, phenotype_name, include_race) {
  frame <- as.data.frame(cd)
  frame$sample_id <- rownames(frame)
  frame$phenotype_raw <- as.numeric(frame[[phenotype_name]])
  frame$age_numeric <- as.numeric(frame$Age)
  frame$pmi_numeric <- as.numeric(frame$PMI)
  frame$n_cells_numeric <- as.numeric(frame$n_cells)
  frame$sex_model <- factor(as.character(frame$Sex), levels = c("F", "M"))
  if (include_race) {
    frame$race_model <- factor(
      as.character(frame$reported_race_group),
      levels = c("White", "Other")
    )
  }
  required <- c(
    "phenotype_raw", "age_numeric", "pmi_numeric", "n_cells_numeric",
    "sex_model"
  )
  if (include_race) required <- c(required, "race_model")
  frame <- droplevels(frame[complete.cases(frame[, required]), , drop = FALSE])
  if (nrow(frame) < 60L) stop("Fewer than 60 complete SEA-AD donors")
  if (length(unique(frame$phenotype_raw)) < 4L) {
    stop("Insufficient pathology levels")
  }
  frame$phenotype_value <- zscore(frame$phenotype_raw)
  frame$age_z <- zscore(frame$age_numeric)
  frame$pmi_z <- zscore(frame$pmi_numeric)
  frame$log10_n_cells_z <- zscore(log10(frame$n_cells_numeric))
  frame
}

fit_one <- function(
    counts, gene_metadata, metadata, phenotype_name, settings, output_dir
) {
  model_formula <- if (settings$include_race) {
    ~ phenotype_value + age_z + sex_model + pmi_z +
      log10_n_cells_z + race_model
  } else {
    ~ phenotype_value + age_z + sex_model + pmi_z + log10_n_cells_z
  }
  design <- model.matrix(model_formula, data = metadata)
  if (qr(design)$rank != ncol(design)) {
    stop("Design is not full rank for ", phenotype_name)
  }
  sample_columns <- match(metadata$sample_id, colnames(counts))
  if (anyNA(sample_columns)) stop("Sample identifiers do not match counts")
  model_counts <- counts[, sample_columns, drop = FALSE]
  storage.mode(model_counts@x) <- "double"
  dge <- DGEList(counts = model_counts)
  minimum_samples <- max(
    settings$min_gene_samples,
    ceiling(settings$min_gene_fraction * ncol(dge))
  )
  keep <- rowSums(cpm(dge) >= settings$cpm_threshold) >= minimum_samples
  if (sum(keep) < 1000L) stop("Implausibly few retained genes")
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  voom_object <- if (settings$quality_weights) {
    voomWithQualityWeights(
      dge, design = design, plot = FALSE, normalize.method = "none"
    )
  } else {
    voom(dge, design = design, plot = FALSE, normalize.method = "none")
  }
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
    gene_metadata[keep, , drop = FALSE],
    data.frame(
      class = settings$cell_class,
      phenotype = phenotype_name,
      beta_D = as.numeric(fit$coefficients[, coefficient]),
      SE_D_moderated = moderated_se,
      t = moderated_t,
      AveExpr = table$AveExpr,
      P.Value = table$P.Value,
      FDR_within_class = table$adj.P.Val,
      B = table$B,
      n_donors = nrow(metadata),
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
  result_path <- file.path(
    output_dir,
    sprintf("%s__%s__gene_results.csv", settings$cell_class, phenotype_name)
  )
  weight_path <- file.path(
    output_dir,
    sprintf("%s__%s__sample_weights.csv", settings$cell_class, phenotype_name)
  )
  write_csv_atomic(result, result_path)
  write_csv_atomic(
    data.frame(
      sample_id = metadata$sample_id,
      donor_id = metadata$donor_id,
      phenotype = phenotype_name,
      phenotype_raw = metadata$phenotype_raw,
      n_cells = metadata$n_cells_numeric,
      sample_quality_weight = sample_weights,
      stringsAsFactors = FALSE
    ),
    weight_path
  )
  data.frame(
    class = settings$cell_class,
    phenotype = phenotype_name,
    status = "complete",
    voom_method = if (settings$quality_weights) {
      "voomWithQualityWeights"
    } else {
      "voom"
    },
    include_race = settings$include_race,
    n_donors = nrow(metadata),
    phenotype_min = min(metadata$phenotype_raw),
    phenotype_max = max(metadata$phenotype_raw),
    genes_tested = nrow(result),
    design_rank = qr(design)$rank,
    design_columns = ncol(design),
    design_condition_number = kappa(design, exact = TRUE),
    minimum_gene_samples = minimum_samples,
    sample_weight_min = min(sample_weights, na.rm = TRUE),
    sample_weight_median = median(sample_weights, na.rm = TRUE),
    sample_weight_max = max(sample_weights, na.rm = TRUE),
    fdr05_genes = sum(result$FDR_within_class < 0.05, na.rm = TRUE),
    result_file = basename(result_path),
    weight_file = basename(weight_path),
    stringsAsFactors = FALSE
  )
}

settings <- parse_args(commandArgs(trailingOnly = TRUE))
input <- normalizePath(settings$input, mustWork = TRUE)
output_dir <- normalizePath(settings$output_dir, mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading ", input)
sce <- readH5AD(input, use_hdf5 = FALSE, reader = "R")
if (!"X" %in% assayNames(sce)) stop("X assay missing")
counts <- assay(sce, "X")
if (!inherits(counts, "sparseMatrix")) counts <- as(counts, "dgCMatrix")
if (nrow(counts) != 36601L) stop("Unexpected gene count")
if (any(counts@x < 0) || any(abs(counts@x - round(counts@x)) > 1e-6)) {
  stop("Counts are not nonnegative integers")
}
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
if (anyDuplicated(gene_metadata$gene_id)) stop("Duplicated gene IDs")

manifest_rows <- list()
for (phenotype_name in settings$phenotypes) {
  result_file <- file.path(
    output_dir,
    sprintf("%s__%s__gene_results.csv", settings$cell_class, phenotype_name)
  )
  if (file.exists(result_file) && !settings$overwrite) {
    message("Skipping existing ", basename(result_file))
    next
  }
  message("Fitting SEA-AD ", settings$cell_class, " / ", phenotype_name)
  metadata <- prepare_metadata(
    colData(sce), phenotype_name, settings$include_race
  )
  manifest_rows[[phenotype_name]] <- fit_one(
    counts, gene_metadata, metadata, phenotype_name, settings, output_dir
  )
}

if (length(manifest_rows)) {
  manifest <- do.call(rbind, manifest_rows)
  manifest_path <- file.path(output_dir, "model_manifest.csv")
  if (file.exists(manifest_path) && !settings$overwrite) {
    old <- read.csv(manifest_path, stringsAsFactors = FALSE)
    manifest <- rbind(old, manifest)
    manifest <- manifest[
      !duplicated(manifest[, c("class", "phenotype")], fromLast = TRUE),
      ,
      drop = FALSE
    ]
  }
  write_csv_atomic(manifest, manifest_path)
}

sink(file.path(output_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()
message("SEA-AD core3 pathology model pipeline complete: ", output_dir)
