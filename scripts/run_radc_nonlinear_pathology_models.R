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
    linear_result_dir = NULL,
    output_dir = NULL,
    min_cells = 20L
  )
  i <- 1L
  while (i <= length(x)) {
    key <- gsub("-", "_", sub("^--", "", x[[i]]))
    if (i == length(x)) stop("Missing value for --", key)
    out[[key]] <- x[[i + 1L]]
    i <- i + 2L
  }
  for (key in c("input_dir", "linear_result_dir", "output_dir")) {
    if (is.null(out[[key]])) stop("--", gsub("_", "-", key), " is required")
  }
  out$min_cells <- as.integer(out$min_cells)
  out
}

zscore <- function(x) {
  value <- as.numeric(scale(as.numeric(x)))
  if (any(!is.finite(value))) stop("Cannot standardize constant covariate")
  value
}

age_numeric <- function(x) {
  as.numeric(sub("\\+$", "", as.character(x)))
}

ancestry_group <- function(x) {
  x <- as.character(x)
  factor(
    ifelse(x %in% c("AFR", "AMR", "EUR"), x, "Other"),
    levels = c("AFR", "AMR", "EUR", "Other")
  )
}

non_ad_diagnosis <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  as.integer(x != "" & x != "AD")
}

pathology_stage <- function(x, phenotype_name) {
  if (phenotype_name == "CERAD") {
    value <- ifelse(x <= 2, "low", ifelse(x == 3, "middle", "high"))
  } else if (phenotype_name == "Braak") {
    value <- ifelse(x <= 2, "low", ifelse(x <= 4, "middle", "high"))
  } else {
    stop("Unknown pathology phenotype")
  }
  factor(value, levels = c("low", "middle", "high"))
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
  frame$phenotype_raw <- if (phenotype_name == "CERAD") {
    as.numeric(frame$CERAD)
  } else if (phenotype_name == "Braak") {
    as.numeric(frame$Braak)
  } else {
    stop("Only CERAD and Braak are supported")
  }
  required <- c(
    "n_cells_numeric", "age_numeric", "pmi_numeric", "sex_model",
    "ancestry_model", "phenotype_raw"
  )
  keep <- complete.cases(frame[, required, drop = FALSE]) &
    frame$n_cells_numeric >= min_cells
  frame <- droplevels(frame[keep, , drop = FALSE])
  frame$age_z <- zscore(frame$age_numeric)
  frame$pmi_z <- zscore(frame$pmi_numeric)
  frame$log10_n_cells_z <- zscore(log10(frame$n_cells_numeric))
  frame$pathology_stage <- pathology_stage(frame$phenotype_raw, phenotype_name)
  polynomial <- poly(frame$phenotype_raw, degree = 2)
  frame$path_linear_basis <- polynomial[, 1]
  frame$path_quadratic_basis <- polynomial[, 2]
  frame
}

moderated_se <- function(fit, coefficient) {
  t_value <- as.numeric(fit$t[, coefficient])
  beta <- as.numeric(fit$coefficients[, coefficient])
  ifelse(is.finite(t_value) & t_value != 0, abs(beta / t_value), NA_real_)
}

settings <- parse_args(commandArgs(trailingOnly = TRUE))
input_dir <- normalizePath(settings$input_dir, mustWork = TRUE)
linear_result_dir <- normalizePath(settings$linear_result_dir, mustWork = TRUE)
output_dir <- normalizePath(settings$output_dir, mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

classes <- c("Astro", "EN", "Endo", "IN", "Immune", "Mural", "OPC", "Oligo")
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
  storage.mode(counts@x) <- "double"
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

  for (phenotype_name in phenotypes) {
    result_path <- file.path(
      output_dir,
      sprintf("%s__%s__nonlinear_results.csv", class_name, phenotype_name)
    )
    if (file.exists(result_path)) {
      message("Skipping existing ", basename(result_path))
      next
    }
    message("Fitting nonlinear ", class_name, " / ", phenotype_name)
    metadata <- prepare_metadata(colData(sce), phenotype_name, settings$min_cells)
    columns <- match(metadata$sample_id, colnames(counts))
    model_counts <- counts[, columns, drop = FALSE]

    polynomial_design <- model.matrix(
      ~ path_linear_basis + path_quadratic_basis + age_z + sex_model + pmi_z +
        ancestry_model + log10_n_cells_z + nonAD_dx,
      data = metadata
    )
    stage_design <- model.matrix(
      ~ pathology_stage + age_z + sex_model + pmi_z + ancestry_model +
        log10_n_cells_z + nonAD_dx,
      data = metadata
    )
    if (
      qr(polynomial_design)$rank != ncol(polynomial_design) ||
      qr(stage_design)$rank != ncol(stage_design)
    ) {
      stop("Nonlinear design is not full rank")
    }

    dge <- DGEList(counts = model_counts)
    minimum_samples <- max(10L, ceiling(0.10 * ncol(dge)))
    keep_gene <- rowSums(cpm(dge) >= 1) >= minimum_samples
    dge <- dge[keep_gene, , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge, method = "TMM")

    voom_poly <- voom(
      dge,
      design = polynomial_design,
      plot = FALSE,
      normalize.method = "none"
    )
    fit_poly <- eBayes(lmFit(voom_poly, polynomial_design), robust = TRUE)
    linear_coef <- which(colnames(polynomial_design) == "path_linear_basis")
    quadratic_coef <- which(colnames(polynomial_design) == "path_quadratic_basis")
    quadratic_table <- topTable(
      fit_poly,
      coef = quadratic_coef,
      number = Inf,
      sort.by = "none"
    )
    omnibus_table <- topTable(
      fit_poly,
      coef = c(linear_coef, quadratic_coef),
      number = Inf,
      sort.by = "none"
    )

    voom_stage <- voom(
      dge,
      design = stage_design,
      plot = FALSE,
      normalize.method = "none"
    )
    fit_stage <- eBayes(lmFit(voom_stage, stage_design), robust = TRUE)
    middle_coef <- which(colnames(stage_design) == "pathology_stagemiddle")
    high_coef <- which(colnames(stage_design) == "pathology_stagehigh")
    middle_table <- topTable(
      fit_stage,
      coef = middle_coef,
      number = Inf,
      sort.by = "none"
    )
    high_table <- topTable(
      fit_stage,
      coef = high_coef,
      number = Inf,
      sort.by = "none"
    )
    stage_omnibus <- topTable(
      fit_stage,
      coef = c(middle_coef, high_coef),
      number = Inf,
      sort.by = "none"
    )

    linear_path <- file.path(
      linear_result_dir,
      sprintf("%s__%s__gene_results.csv", class_name, phenotype_name)
    )
    linear_primary <- read.csv(linear_path, stringsAsFactors = FALSE)[
      , c("gene_id", "beta_D", "SE_D_moderated", "P.Value", "FDR_within_class")
    ]
    linear_primary <- linear_primary[match(
      gene_metadata$gene_id[keep_gene],
      linear_primary$gene_id
    ), ]
    if (anyNA(linear_primary$gene_id)) stop("Linear-primary gene match failed")

    result <- cbind(
      gene_metadata[keep_gene, , drop = FALSE],
      data.frame(
        class = class_name,
        phenotype = phenotype_name,
        n_donors = nrow(metadata),
        primary_linear_beta_D = linear_primary$beta_D,
        primary_linear_SE_D = linear_primary$SE_D_moderated,
        primary_linear_P = linear_primary$P.Value,
        primary_linear_FDR = linear_primary$FDR_within_class,
        quadratic_basis_beta = as.numeric(
          fit_poly$coefficients[, quadratic_coef]
        ),
        quadratic_basis_SE = moderated_se(fit_poly, quadratic_coef),
        quadratic_P = quadratic_table$P.Value,
        quadratic_FDR = quadratic_table$adj.P.Val,
        polynomial_omnibus_F = omnibus_table$F,
        polynomial_omnibus_P = omnibus_table$P.Value,
        polynomial_omnibus_FDR = omnibus_table$adj.P.Val,
        middle_vs_low_beta = as.numeric(fit_stage$coefficients[, middle_coef]),
        middle_vs_low_SE = moderated_se(fit_stage, middle_coef),
        middle_vs_low_P = middle_table$P.Value,
        middle_vs_low_FDR = middle_table$adj.P.Val,
        high_vs_low_beta = as.numeric(fit_stage$coefficients[, high_coef]),
        high_vs_low_SE = moderated_se(fit_stage, high_coef),
        high_vs_low_P = high_table$P.Value,
        high_vs_low_FDR = high_table$adj.P.Val,
        stage_omnibus_F = stage_omnibus$F,
        stage_omnibus_P = stage_omnibus$P.Value,
        stage_omnibus_FDR = stage_omnibus$adj.P.Val,
        stringsAsFactors = FALSE
      )
    )
    write.csv(result, result_path, row.names = FALSE)
    stage_counts <- table(metadata$pathology_stage)
    manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
      class = class_name,
      phenotype = phenotype_name,
      n_donors = nrow(metadata),
      low_n = unname(stage_counts[["low"]]),
      middle_n = unname(stage_counts[["middle"]]),
      high_n = unname(stage_counts[["high"]]),
      genes = nrow(result),
      quadratic_fdr05 = sum(result$quadratic_FDR < 0.05),
      stage_omnibus_fdr05 = sum(result$stage_omnibus_FDR < 0.05),
      high_vs_low_fdr05 = sum(result$high_vs_low_FDR < 0.05),
      linear_high_low_beta_spearman = cor(
        result$primary_linear_beta_D,
        result$high_vs_low_beta,
        method = "spearman"
      ),
      prioritized_sign_concordance = {
        selected <- result$primary_linear_P < 0.01 |
          result$high_vs_low_P < 0.01
        if (any(selected)) {
          mean(
            sign(result$primary_linear_beta_D[selected]) ==
              sign(result$high_vs_low_beta[selected])
          )
        } else {
          NA_real_
        }
      },
      result_file = basename(result_path),
      stringsAsFactors = FALSE
    )
    rm(
      model_counts, dge, voom_poly, fit_poly, voom_stage, fit_stage, result
    )
    gc()
  }
  rm(sce, counts)
  gc()
}

manifest <- do.call(rbind, manifest_rows)
write.csv(manifest, file.path(output_dir, "nonlinear_model_manifest.csv"), row.names = FALSE)
sink(file.path(output_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()
message("RADC nonlinear pathology models complete: ", output_dir)
