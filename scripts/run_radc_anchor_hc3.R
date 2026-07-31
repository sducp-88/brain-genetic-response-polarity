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
    primary_dir = NULL,
    anchors = NULL,
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
  for (key in c("input_dir", "primary_dir", "anchors", "output_dir")) {
    if (is.null(out[[key]])) stop("--", gsub("_", "-", key), " is required")
  }
  out$min_cells <- as.integer(out$min_cells)
  out
}

clean_gene_id <- function(x) sub("\\..*$", "", as.character(x))

zscore <- function(x) {
  value <- as.numeric(scale(as.numeric(x)))
  if (any(!is.finite(value))) stop("Cannot standardize constant covariate")
  value
}

age_numeric <- function(x) as.numeric(sub("\\+$", "", as.character(x)))

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
  } else {
    as.numeric(frame$Braak)
  }
  required <- c(
    "n_cells_numeric", "age_numeric", "pmi_numeric", "sex_model",
    "ancestry_model", "phenotype_raw"
  )
  keep <- complete.cases(frame[, required, drop = FALSE]) &
    frame$n_cells_numeric >= min_cells
  frame <- droplevels(frame[keep, , drop = FALSE])
  frame$phenotype_value <- zscore(frame$phenotype_raw)
  frame$age_z <- zscore(frame$age_numeric)
  frame$pmi_z <- zscore(frame$pmi_numeric)
  frame$log10_n_cells_z <- zscore(log10(frame$n_cells_numeric))
  frame
}

hc3_for_gene <- function(y, weights, design) {
  if (any(!is.finite(y)) || any(!is.finite(weights)) || any(weights <= 0)) {
    return(c(beta = NA_real_, se = NA_real_, max_hat = NA_real_))
  }
  sqrt_w <- sqrt(weights)
  xw <- design * sqrt_w
  yw <- y * sqrt_w
  inverse <- tryCatch(
    solve(crossprod(xw)),
    error = function(e) NULL
  )
  if (is.null(inverse)) {
    return(c(beta = NA_real_, se = NA_real_, max_hat = NA_real_))
  }
  beta <- as.numeric(inverse %*% crossprod(xw, yw))
  residual <- y - as.numeric(design %*% beta)
  leverage <- rowSums((xw %*% inverse) * xw)
  if (any(leverage >= 1 - 1e-8)) {
    return(c(beta = NA_real_, se = NA_real_, max_hat = max(leverage)))
  }
  omega <- (sqrt_w * residual / (1 - leverage))^2
  meat <- crossprod(xw, xw * omega)
  covariance <- inverse %*% meat %*% inverse
  coefficient <- which(colnames(design) == "phenotype_value")
  c(
    beta = beta[[coefficient]],
    se = sqrt(covariance[coefficient, coefficient]),
    max_hat = max(leverage)
  )
}

write_csv_atomic <- function(x, path) {
  temporary <- paste0(path, ".tmp")
  write.csv(x, temporary, row.names = FALSE, quote = TRUE, na = "")
  if (file.exists(path)) file.remove(path)
  if (!file.rename(temporary, path)) stop("Atomic rename failed")
}

settings <- parse_args(commandArgs(trailingOnly = TRUE))
input_dir <- normalizePath(settings$input_dir, mustWork = TRUE)
primary_dir <- normalizePath(settings$primary_dir, mustWork = TRUE)
output_dir <- normalizePath(settings$output_dir, mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

anchors <- read.csv(settings$anchors, stringsAsFactors = FALSE)
anchors <- anchors[
  anchors$disease == "AD" &
    anchors$evidence_grade %in% c("G1", "G2") &
    anchors$primary_anchor_eligible == "yes",
  ,
  drop = FALSE
]
anchors$gene_id_clean <- clean_gene_id(anchors$gene_id)
classes <- sort(unique(anchors$cell_class))
phenotypes <- c("CERAD", "Braak")
audit_rows <- list()

for (class_name in classes) {
  message("Reading ", class_name)
  input_path <- file.path(
    input_dir, sprintf("RADC_%s_pseudobulk_counts.h5ad", class_name)
  )
  sce <- readH5AD(input_path, use_hdf5 = FALSE, reader = "R")
  counts <- assay(sce, "X")
  if (!inherits(counts, "sparseMatrix")) counts <- as(counts, "dgCMatrix")

  for (phenotype_name in phenotypes) {
    message("HC3 ", class_name, " / ", phenotype_name)
    metadata <- prepare_metadata(colData(sce), phenotype_name, settings$min_cells)
    design <- model.matrix(
      ~ phenotype_value + age_z + sex_model + pmi_z + ancestry_model +
        log10_n_cells_z + nonAD_dx,
      data = metadata
    )
    if (qr(design)$rank != ncol(design)) stop("Design is not full rank")
    columns <- match(metadata$sample_id, colnames(counts))
    model_counts <- counts[, columns, drop = FALSE]
    storage.mode(model_counts@x) <- "double"
    dge <- DGEList(counts = model_counts)
    minimum_samples <- max(10L, ceiling(0.10 * ncol(dge)))
    keep_gene <- rowSums(cpm(dge) >= 1) >= minimum_samples
    dge <- dge[keep_gene, , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge, method = "TMM")
    voom_object <- voom(
      dge, design = design, plot = FALSE, normalize.method = "none"
    )
    clean_ids <- clean_gene_id(rownames(voom_object$E))
    class_anchor_ids <- unique(
      anchors$gene_id_clean[anchors$cell_class == class_name]
    )
    rows <- match(class_anchor_ids, clean_ids)
    estimable <- !is.na(rows)
    primary_path <- file.path(
      primary_dir,
      sprintf("%s__%s__gene_results.csv", class_name, phenotype_name)
    )
    primary <- read.csv(primary_path, stringsAsFactors = FALSE)
    primary$gene_id_clean <- clean_gene_id(primary$gene_id)

    result_rows <- lapply(seq_along(class_anchor_ids), function(index) {
      gene_id <- class_anchor_ids[[index]]
      if (!estimable[[index]]) {
        return(data.frame(
          gene_id_clean = gene_id,
          beta_D = NA_real_,
          SE_D_moderated = NA_real_,
          AveExpr = NA_real_,
          beta_D_HC3 = NA_real_,
          SE_D_HC3 = NA_real_,
          max_hat = NA_real_
        ))
      }
      row <- rows[[index]]
      hc3 <- hc3_for_gene(
        as.numeric(voom_object$E[row, ]),
        as.numeric(voom_object$weights[row, ]),
        design
      )
      primary_row <- primary[primary$gene_id_clean == gene_id, , drop = FALSE]
      if (nrow(primary_row) != 1L) stop("Primary gene match failed")
      if (
        is.finite(hc3[["beta"]]) &&
        abs(hc3[["beta"]] - primary_row$beta_D) > 1e-6
      ) {
        stop("HC3 WLS coefficient does not reproduce primary coefficient")
      }
      data.frame(
        gene_id_clean = gene_id,
        beta_D = primary_row$beta_D,
        SE_D_moderated = primary_row$SE_D_moderated,
        AveExpr = primary_row$AveExpr,
        beta_D_HC3 = hc3[["beta"]],
        SE_D_HC3 = hc3[["se"]],
        max_hat = hc3[["max_hat"]]
      )
    })
    result <- do.call(rbind, result_rows)
    result$class <- class_name
    result$phenotype <- phenotype_name
    result$n_donors <- nrow(metadata)
    output_path <- file.path(
      output_dir, sprintf("%s__%s__anchor_hc3.csv", class_name, phenotype_name)
    )
    write_csv_atomic(result, output_path)
    audit_rows[[length(audit_rows) + 1L]] <- data.frame(
      class = class_name,
      phenotype = phenotype_name,
      n_donors = nrow(metadata),
      anchor_genes = nrow(result),
      hc3_estimable = sum(
        is.finite(result$beta_D_HC3) & is.finite(result$SE_D_HC3) &
          result$SE_D_HC3 > 0
      ),
      max_hat = max(result$max_hat, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  rm(sce, counts)
  gc()
}

write_csv_atomic(
  do.call(rbind, audit_rows), file.path(output_dir, "hc3_model_manifest.csv")
)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("COMPLETE RADC anchor HC3")

