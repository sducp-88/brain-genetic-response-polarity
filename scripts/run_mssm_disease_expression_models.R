#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(SingleCellExperiment)
  library(splines)
  library(zellkonverter)
})

parse_arguments <- function(arguments) {
  defaults <- list(
    pseudobulk_dir = NA_character_,
    qc_dir = NA_character_,
    readiness = NA_character_,
    weights_dir = NA_character_,
    anchors = NA_character_,
    output_dir = NA_character_
  )
  index <- 1L
  while (index <= length(arguments)) {
    key <- sub("^--", "", arguments[[index]])
    if (!key %in% names(defaults) || index == length(arguments)) {
      stop("Unknown or incomplete argument: ", arguments[[index]])
    }
    defaults[[key]] <- arguments[[index + 1L]]
    index <- index + 2L
  }
  missing <- names(defaults)[is.na(unlist(defaults))]
  if (length(missing)) {
    stop("Missing required arguments: ", paste(missing, collapse = ", "))
  }
  defaults
}

sha256_file <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("The digest package is required for SHA-256 manifests.")
  }
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

write_json_atomic <- function(payload, path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("The jsonlite package is required.")
  }
  temporary <- paste0(path, ".tmp")
  writeLines(
    jsonlite::toJSON(
      payload,
      auto_unbox = TRUE,
      pretty = TRUE,
      na = "null",
      null = "null"
    ),
    temporary,
    useBytes = TRUE
  )
  if (!file.rename(temporary, path)) {
    stop("Could not atomically write JSON: ", path)
  }
}

clean_gene_id <- function(values) {
  sub("\\..*$", "", as.character(values))
}

required_columns <- function(frame, columns, label) {
  missing <- setdiff(columns, names(frame))
  if (length(missing)) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  }
}

effective_sample_size <- function(weights) {
  sum(weights)^2 / sum(weights^2)
}

weighted_hc3 <- function(response, design, weights, coefficient_index) {
  weighted_fit <- lm.wfit(
    x = design,
    y = response,
    w = weights
  )
  bread <- solve(crossprod(design, design * weights))
  leverage <- weights * rowSums((design %*% bread) * design)
  leverage <- pmin(leverage, 1 - 1e-8)
  score_scale <- (
    weights * weighted_fit$residuals / pmax(1 - leverage, 1e-8)
  )
  meat <- crossprod(design, design * score_scale^2)
  covariance <- bread %*% meat %*% bread
  beta <- weighted_fit$coefficients[[coefficient_index]]
  standard_error <- sqrt(covariance[coefficient_index, coefficient_index])
  c(
    beta = beta,
    SE = standard_error,
    p_value = 2 * pnorm(-abs(beta / standard_error))
  )
}

arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))
dir.create(arguments$output_dir, recursive = TRUE, showWarnings = FALSE)
result_directory <- file.path(arguments$output_dir, "all_gene_results")
dir.create(result_directory, recursive = TRUE, showWarnings = FALSE)

readiness <- fread(arguments$readiness)
required_columns(
  readiness,
  c(
    "disease",
    "class",
    "cases",
    "controls",
    "confirmatory_ready",
    "case_ess",
    "control_ess"
  ),
  "Readiness table"
)
readiness[, confirmatory_ready := as.logical(confirmatory_ready)]
ready_rows <- readiness[confirmatory_ready == TRUE]
if (!nrow(ready_rows)) {
  stop("No disease-by-cell-class comparison passed readiness gates.")
}

anchors <- fread(arguments$anchors)
required_columns(
  anchors,
  c(
    "anchor_unit_id",
    "disease",
    "gene_id",
    "gene_symbol",
    "cell_class",
    "beta_G",
    "SE_G",
    "evidence_grade",
    "primary_anchor_eligible"
  ),
  "Genetic anchor table"
)
anchors <- anchors[
  primary_anchor_eligible == "yes" &
    evidence_grade %in% c("G1", "G2")
]
anchors[, gene_id_clean := clean_gene_id(gene_id)]

diagnostic_rows <- list()
all_result_rows <- list()
anchor_result_rows <- list()
output_hashes <- list()

for (row_index in seq_len(nrow(ready_rows))) {
  disease_name <- ready_rows$disease[[row_index]]
  class_name <- ready_rows$class[[row_index]]
  comparison_label <- paste(disease_name, class_name, sep = "__")
  message("Starting ", comparison_label)

  h5ad_path <- file.path(
    arguments$pseudobulk_dir,
    paste0("MSSM_", class_name, "_pseudobulk_counts.h5ad")
  )
  filter_path <- file.path(
    arguments$qc_dir,
    paste0(class_name, "_gene_filter.csv")
  )
  weight_path <- file.path(
    arguments$weights_dir,
    paste0(disease_name, "_", class_name, "_overlap_weights.csv")
  )
  if (!file.exists(h5ad_path) || !file.exists(filter_path) ||
      !file.exists(weight_path)) {
    stop("Required input is absent for ", comparison_label)
  }

  weights <- fread(weight_path)
  required_columns(
    weights,
    c(
      "donor_id",
      "analysis_role",
      "age_numeric",
      "Sex",
      "ancestry_harmonized",
      "PMI_numeric",
      "overlap_weight"
    ),
    paste0("Weights for ", comparison_label)
  )
  if (anyDuplicated(weights$donor_id)) {
    stop("Duplicate donor in weight file: ", comparison_label)
  }
  weights[, case := as.integer(analysis_role == "case")]
  weights[, Sex := factor(Sex)]
  weights[, ancestry_harmonized := factor(ancestry_harmonized)]
  weights[, overlap_weight_normalized := overlap_weight / mean(overlap_weight)]

  sce <- readH5AD(h5ad_path, use_hdf5 = FALSE, reader = "R")
  assay_name <- if ("X" %in% assayNames(sce)) "X" else assayNames(sce)[[1]]
  donor_ids <- as.character(colData(sce)$donor_id)
  if (anyDuplicated(donor_ids)) {
    stop("Duplicate donor columns in ", h5ad_path)
  }
  sample_indices <- match(weights$donor_id, donor_ids)
  if (anyNA(sample_indices)) {
    stop("Weight donors absent from H5AD: ", comparison_label)
  }
  sce <- sce[, sample_indices]
  if (!identical(as.character(colData(sce)$donor_id), weights$donor_id)) {
    stop("Donor ordering failed: ", comparison_label)
  }

  gene_ids <- rownames(sce)
  gene_ids_clean <- clean_gene_id(gene_ids)
  if (anyDuplicated(gene_ids_clean)) {
    stop("Cleaned gene identifiers are not unique in ", h5ad_path)
  }
  gene_metadata <- as.data.table(as.data.frame(rowData(sce)))
  gene_metadata[, gene_id := gene_ids]
  gene_metadata[, gene_id_clean := gene_ids_clean]
  if (!"gene_name" %in% names(gene_metadata)) {
    gene_metadata[, gene_name := NA_character_]
  }

  gene_filter <- fread(filter_path)
  required_columns(
    gene_filter,
    c("gene_id", "retain_primary_filter"),
    paste0("Gene filter for ", class_name)
  )
  gene_filter[, gene_id_clean := clean_gene_id(gene_id)]
  retained_gene_ids <- gene_filter[
    as.logical(retain_primary_filter) == TRUE,
    unique(gene_id_clean)
  ]
  keep <- gene_ids_clean %in% retained_gene_ids
  if (sum(keep) < 1000L) {
    stop("Too few genes pass the frozen filter: ", comparison_label)
  }

  counts <- as.matrix(assay(sce, assay_name)[keep, , drop = FALSE])
  storage.mode(counts) <- "integer"
  if (any(!is.finite(counts)) || any(counts < 0L)) {
    stop("Invalid counts in ", comparison_label)
  }
  rownames(counts) <- gene_ids_clean[keep]
  colnames(counts) <- weights$donor_id

  design <- model.matrix(
    ~ case +
      ns(age_numeric, df = 3) +
      Sex +
      ancestry_harmonized +
      PMI_numeric,
    data = weights
  )
  if (qr(design)$rank != ncol(design)) {
    stop("Rank-deficient design for ", comparison_label)
  }
  rownames(design) <- weights$donor_id
  if (!"case" %in% colnames(design)) {
    stop("Disease coefficient was not constructed: ", comparison_label)
  }

  dge <- DGEList(counts = counts)
  dge <- calcNormFactors(dge, method = "TMM")
  voom_object <- voomWithQualityWeights(
    dge,
    design = design,
    plot = FALSE,
    normalize.method = "none"
  )
  voom_object$weights <- sweep(
    voom_object$weights,
    MARGIN = 2L,
    STATS = weights$overlap_weight_normalized,
    FUN = "*"
  )
  fit <- lmFit(voom_object, design)
  fit <- eBayes(fit, robust = TRUE)
  coefficient_index <- match("case", colnames(fit$coefficients))
  beta_d <- fit$coefficients[, coefficient_index]
  se_d <- fit$stdev.unscaled[, coefficient_index] * sqrt(fit$s2.post)
  p_value <- fit$p.value[, coefficient_index]
  result <- data.table(
    disease = disease_name,
    class = class_name,
    gene_id_clean = rownames(fit$coefficients),
    beta_D = as.numeric(beta_d),
    SE_D = as.numeric(se_d),
    t_moderated = as.numeric(fit$t[, coefficient_index]),
    p_value = as.numeric(p_value),
    FDR = p.adjust(p_value, method = "BH"),
    average_log2_expression = as.numeric(fit$Amean),
    cases = sum(weights$case == 1L),
    controls = sum(weights$case == 0L),
    case_ESS = effective_sample_size(
      weights[case == 1L, overlap_weight]
    ),
    control_ESS = effective_sample_size(
      weights[case == 0L, overlap_weight]
    ),
    model = paste(
      "voomWithQualityWeights + normalized overlap weights;",
      "case + ns(age,3) + sex + harmonized ancestry + PMI"
    )
  )
  result[
    ,
    `:=`(
      beta_D_HC3 = NA_real_,
      SE_D_HC3 = NA_real_,
      p_value_HC3 = NA_real_
    )
  ]
  anchor_gene_ids <- unique(
    anchors[
      disease == disease_name & cell_class == class_name,
      gene_id_clean
    ]
  )
  anchor_gene_ids <- intersect(
    anchor_gene_ids,
    rownames(voom_object$E)
  )
  for (anchor_gene_id in anchor_gene_ids) {
    anchor_gene_index <- match(
      anchor_gene_id,
      rownames(voom_object$E)
    )
    if (is.na(anchor_gene_index)) {
      next
    }
    robust_result <- weighted_hc3(
      response = voom_object$E[anchor_gene_index, ],
      design = design,
      weights = voom_object$weights[anchor_gene_index, ],
      coefficient_index = coefficient_index
    )
    result[
      gene_id_clean == anchor_gene_id,
      `:=`(
        beta_D_HC3 = robust_result[["beta"]],
        SE_D_HC3 = robust_result[["SE"]],
        p_value_HC3 = robust_result[["p_value"]]
      )
    ]
  }
  result <- merge(
    result,
    gene_metadata[, .(gene_id, gene_id_clean, gene_name)],
    by = "gene_id_clean",
    all.x = TRUE,
    sort = FALSE
  )
  setcolorder(
    result,
    c(
      "disease",
      "class",
      "gene_id",
      "gene_id_clean",
      "gene_name",
      setdiff(names(result), c(
        "disease",
        "class",
        "gene_id",
        "gene_id_clean",
        "gene_name"
      ))
    )
  )
  result_path <- file.path(
    result_directory,
    paste0(comparison_label, "_all_gene_results.csv.gz")
  )
  fwrite(result, result_path)
  output_hashes[[basename(result_path)]] <- sha256_file(result_path)
  all_result_rows[[comparison_label]] <- result

  comparison_anchors <- anchors[
    disease == disease_name & cell_class == class_name
  ]
  if (nrow(comparison_anchors)) {
    result_for_anchor <- result[
      ,
      setdiff(names(result), c("disease", "gene_id")),
      with = FALSE
    ]
    anchor_output <- merge(
      comparison_anchors,
      result_for_anchor,
      by = "gene_id_clean",
      all.x = TRUE,
      sort = FALSE,
      allow.cartesian = TRUE
    )
    anchor_output[, expression_estimable := !is.na(beta_D)]
    anchor_result_rows[[comparison_label]] <- anchor_output
  }

  diagnostic_rows[[comparison_label]] <- data.table(
    disease = disease_name,
    class = class_name,
    h5ad_sha256 = sha256_file(h5ad_path),
    weight_sha256 = sha256_file(weight_path),
    filter_sha256 = sha256_file(filter_path),
    samples = nrow(weights),
    cases = sum(weights$case == 1L),
    controls = sum(weights$case == 0L),
    case_ESS = effective_sample_size(
      weights[case == 1L, overlap_weight]
    ),
    control_ESS = effective_sample_size(
      weights[case == 0L, overlap_weight]
    ),
    design_columns = ncol(design),
    design_rank = qr(design)$rank,
    genes_retained = nrow(result),
    significant_genes_FDR_0_05 = sum(result$FDR < 0.05, na.rm = TRUE)
  )

  rm(
    sce,
    counts,
    dge,
    voom_object,
    fit,
    result
  )
  invisible(gc())
  message("Completed ", comparison_label)
}

diagnostics <- rbindlist(diagnostic_rows, use.names = TRUE, fill = TRUE)
diagnostics_path <- file.path(
  arguments$output_dir,
  "mssm_disease_expression_model_diagnostics.csv"
)
fwrite(diagnostics, diagnostics_path)
output_hashes[[basename(diagnostics_path)]] <- sha256_file(diagnostics_path)

anchor_results <- if (length(anchor_result_rows)) {
  rbindlist(anchor_result_rows, use.names = TRUE, fill = TRUE)
} else {
  copy(anchors)[0]
}
if (nrow(anchor_results)) {
  anchor_results[
    ,
    SE_D_direction := fifelse(
      is.finite(SE_D_HC3) & SE_D_HC3 > 0,
      SE_D_HC3,
      SE_D
    )
  ]
  anchor_results[
    !is.na(beta_D) & is.finite(SE_D_direction) & SE_D_direction > 0 &
      is.finite(beta_G) & is.finite(SE_G) & SE_G > 0,
    `:=`(
      P_beta_D_positive = pnorm(beta_D / SE_D_direction),
      P_beta_G_positive_model = pnorm(beta_G / SE_G)
    )
  ]
  anchor_results[
    !is.na(P_beta_D_positive) & !is.na(P_beta_G_positive_model),
    P_aligned_independence := (
      P_beta_D_positive * P_beta_G_positive_model +
        (1 - P_beta_D_positive) * (1 - P_beta_G_positive_model)
    )
  ]
  anchor_results[
    !is.na(P_aligned_independence),
    S_independence := 2 * P_aligned_independence - 1
  ]
}
anchor_path <- file.path(
  arguments$output_dir,
  "mssm_G1_G2_anchor_disease_effects.csv"
)
fwrite(anchor_results, anchor_path)
output_hashes[[basename(anchor_path)]] <- sha256_file(anchor_path)

manifest <- list(
  status = "COMPLETE",
  expression_effects_estimated = TRUE,
  primary_model = paste(
    "voomWithQualityWeights; normalized overlap weights;",
    "case + ns(age,3) + sex + harmonized ancestry + PMI"
  ),
  anchor_uncertainty = paste(
    "HC3 sandwich SE on the same weighted log-CPM model is primary for",
    "directional coupling; moderated limma SE is retained for comparison"
  ),
  gene_filter = paste(
    "Frozen QC filter: CPM >= 1 in >= max(10, 10% of eligible",
    "pseudobulk samples)"
  ),
  disease_cell_comparisons = nrow(diagnostics),
  anchor_rows = nrow(anchor_results),
  estimable_anchor_rows = if (nrow(anchor_results)) {
    sum(anchor_results$expression_estimable, na.rm = TRUE)
  } else {
    0L
  },
  caution = paste(
    "P_aligned_independence assumes independent genetic and disease-effect",
    "errors; rho-grid sensitivity is required before strongest inference."
  ),
  output_sha256 = output_hashes,
  completed_at = format(
    Sys.time(),
    "%Y-%m-%dT%H:%M:%S%z",
    tz = "Asia/Shanghai"
  )
)
write_json_atomic(
  manifest,
  file.path(
    arguments$output_dir,
    "mssm_disease_expression_model_manifest.json"
  )
)
cat(
  "Completed",
  nrow(diagnostics),
  "disease-by-cell-class models and",
  nrow(anchor_results),
  "anchor rows.\n"
)
