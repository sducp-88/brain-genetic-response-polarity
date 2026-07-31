#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(zellkonverter)
  library(edgeR)
  library(limma)
  library(Matrix)
})

FROZEN_REGIONS <- c("DFC", "MEC", "MTG")

parse_args <- function(x) {
  out <- list(
    input = NULL,
    output_dir = NULL,
    cell_class = NULL,
    analysis = "region",
    region = "DFC",
    phenotypes = "CERAD,Braak",
    composition_pcs = FALSE,
    cpm_threshold = 1,
    min_gene_samples = 10L,
    min_gene_fraction = 0.10,
    overwrite = FALSE
  )
  i <- 1L
  while (i <= length(x)) {
    key <- sub("^--", "", x[[i]])
    if (key %in% c("composition-pcs", "overwrite")) {
      out[[gsub("-", "_", key)]] <- TRUE
      i <- i + 1L
    } else {
      if (i == length(x)) stop("Missing value for --", key)
      out[[gsub("-", "_", key)]] <- x[[i + 1L]]
      i <- i + 2L
    }
  }
  if (
    is.null(out$input) ||
    is.null(out$output_dir) ||
    is.null(out$cell_class)
  ) {
    stop("--input, --output-dir, and --cell-class are required")
  }
  if (!out$analysis %in% c("region", "repeated")) {
    stop("--analysis must be region or repeated")
  }
  if (!out$region %in% FROZEN_REGIONS) {
    stop("--region must be DFC, MEC, or MTG")
  }
  out$phenotypes <- trimws(
    strsplit(out$phenotypes, ",", fixed = TRUE)[[1]]
  )
  if (
    !length(out$phenotypes) ||
    any(!out$phenotypes %in% c("CERAD", "Braak")) ||
    anyDuplicated(out$phenotypes)
  ) {
    stop("--phenotypes must be a unique subset of CERAD,Braak")
  }
  out$cpm_threshold <- as.numeric(out$cpm_threshold)
  out$min_gene_samples <- as.integer(out$min_gene_samples)
  out$min_gene_fraction <- as.numeric(out$min_gene_fraction)
  out
}

write_csv_atomic <- function(x, path) {
  temporary <- paste0(path, ".tmp")
  write.csv(x, temporary, row.names = FALSE, quote = TRUE, na = "")
  if (file.exists(path)) file.remove(path)
  if (!file.rename(temporary, path)) {
    stop("Failed atomic rename: ", path)
  }
}

zscore <- function(x, label) {
  result <- as.numeric(scale(as.numeric(x)))
  if (any(!is.finite(result))) {
    stop("Cannot standardize ", label)
  }
  result
}

donor_zscore <- function(frame, column) {
  unique_frame <- frame[
    !duplicated(frame$donor_id),
    c("donor_id", column),
    drop = FALSE
  ]
  standardized <- zscore(unique_frame[[column]], column)
  lookup <- setNames(standardized, unique_frame$donor_id)
  as.numeric(lookup[frame$donor_id])
}

prepare_metadata <- function(cd, phenotype_name, settings) {
  frame <- as.data.frame(cd)
  frame$sample_id <- rownames(frame)
  frame$donor_id <- as.character(frame$donor_id)
  frame$region_model <- factor(
    as.character(frame$region),
    levels = FROZEN_REGIONS
  )
  if (settings$analysis == "region") {
    frame <- frame[
      as.character(frame$region_model) == settings$region,
      ,
      drop = FALSE
    ]
    frame$region_model <- droplevels(frame$region_model)
  }
  frame$phenotype_raw <- as.numeric(frame[[phenotype_name]])
  frame$age_numeric <- as.numeric(frame$Age)
  frame$pmi_numeric <- as.numeric(frame$PMI)
  frame$n_cells_numeric <- as.numeric(frame$n_cells)
  frame$sex_model <- factor(as.character(frame$Sex), levels = c("F", "M"))

  pc_columns <- grep(
    "^composition_PC[1-3]$",
    colnames(frame),
    value = TRUE
  )
  pc_columns <- pc_columns[order(pc_columns)]
  if (settings$composition_pcs && !length(pc_columns)) {
    stop("Composition-PC model requested but no PCs are present")
  }
  required <- c(
    "phenotype_raw",
    "age_numeric",
    "pmi_numeric",
    "n_cells_numeric",
    "sex_model",
    "donor_id",
    "region_model"
  )
  if (settings$composition_pcs) {
    required <- c(required, pc_columns)
  }
  frame <- droplevels(
    frame[complete.cases(frame[, required, drop = FALSE]), , drop = FALSE]
  )

  if (settings$analysis == "repeated") {
    region_count <- tapply(
      as.character(frame$region_model),
      frame$donor_id,
      function(value) length(unique(value))
    )
    retained_donors <- names(region_count)[region_count >= 2L]
    frame <- droplevels(
      frame[frame$donor_id %in% retained_donors, , drop = FALSE]
    )
  }
  if (settings$composition_pcs) {
    variable_pc <- vapply(
      frame[, pc_columns, drop = FALSE],
      function(value) {
        observed_sd <- sd(as.numeric(value), na.rm = TRUE)
        is.finite(observed_sd) && observed_sd > 1e-10
      },
      logical(1)
    )
    pc_columns <- pc_columns[variable_pc]
    if (!length(pc_columns)) {
      stop("No variable composition PC remains in the analysis subset")
    }
  }
  if (nrow(frame) < 60L || length(unique(frame$donor_id)) < 60L) {
    stop("Fewer than 60 complete SEA-AD donors")
  }
  if (length(unique(frame$phenotype_raw)) < 4L) {
    stop("Insufficient pathology levels")
  }

  if (settings$analysis == "repeated") {
    frame$phenotype_value <- donor_zscore(frame, "phenotype_raw")
    frame$age_z <- donor_zscore(frame, "age_numeric")
    frame$pmi_z <- donor_zscore(frame, "pmi_numeric")
  } else {
    frame$phenotype_value <- zscore(
      frame$phenotype_raw, "phenotype_raw"
    )
    frame$age_z <- zscore(frame$age_numeric, "age_numeric")
    frame$pmi_z <- zscore(frame$pmi_numeric, "pmi_numeric")
  }
  frame$log10_n_cells_z <- zscore(
    log10(frame$n_cells_numeric),
    "log10_n_cells"
  )
  if (settings$composition_pcs) {
    for (column in pc_columns) {
      frame[[paste0(column, "_z")]] <- zscore(frame[[column]], column)
    }
  }
  attr(frame, "composition_columns") <- if (settings$composition_pcs) {
    paste0(pc_columns, "_z")
  } else {
    character()
  }
  frame
}

build_design <- function(metadata, settings) {
  composition_columns <- attr(metadata, "composition_columns")
  terms <- c(
    if (settings$analysis == "repeated") {
      "phenotype_value * region_model"
    } else {
      "phenotype_value"
    },
    "age_z",
    "sex_model",
    "pmi_z",
    "log10_n_cells_z",
    composition_columns
  )
  formula <- as.formula(
    paste("~", paste(terms, collapse = " + "))
  )
  design <- model.matrix(formula, data = metadata)
  if (qr(design)$rank != ncol(design)) {
    stop(
      "Design is not full rank: rank ",
      qr(design)$rank,
      " vs ",
      ncol(design)
    )
  }
  list(
    design = design,
    formula = paste(deparse(formula), collapse = ""),
    composition_columns = composition_columns
  )
}

fit_one <- function(
    counts,
    gene_metadata,
    metadata,
    phenotype_name,
    settings,
    output_dir
) {
  design_record <- build_design(metadata, settings)
  design <- design_record$design
  sample_columns <- match(metadata$sample_id, colnames(counts))
  if (anyNA(sample_columns)) {
    stop("Sample identifiers do not match counts")
  }
  model_counts <- counts[, sample_columns, drop = FALSE]
  storage.mode(model_counts@x) <- "double"
  dge <- DGEList(counts = model_counts)
  minimum_samples <- max(
    settings$min_gene_samples,
    ceiling(settings$min_gene_fraction * ncol(dge))
  )
  keep <- rowSums(cpm(dge) >= settings$cpm_threshold) >= minimum_samples
  if (sum(keep) < 1000L) {
    stop("Implausibly few retained genes")
  }
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  voom_object <- voom(
    dge,
    design = design,
    plot = FALSE,
    normalize.method = "none"
  )

  consensus_correlation <- NA_real_
  if (settings$analysis == "repeated") {
    correlation_fit <- duplicateCorrelation(
      voom_object,
      design,
      block = metadata$donor_id
    )
    consensus_correlation <- correlation_fit$consensus.correlation
    if (!is.finite(consensus_correlation)) {
      stop("Non-finite duplicateCorrelation estimate")
    }
    base_fit <- lmFit(
      voom_object,
      design,
      block = metadata$donor_id,
      correlation = consensus_correlation
    )
  } else {
    base_fit <- lmFit(voom_object, design)
  }
  region_contrast_table <- NULL
  interaction_table <- NULL
  main_effect_label <- "region_specific_pathology"
  if (settings$analysis == "repeated") {
    coefficient_names <- colnames(design)
    phenotype_index <- which(coefficient_names == "phenotype_value")
    mec_interaction_index <- which(
      coefficient_names %in% c(
        "phenotype_value:region_modelMEC",
        "region_modelMEC:phenotype_value"
      )
    )
    mtg_interaction_index <- which(
      coefficient_names %in% c(
        "phenotype_value:region_modelMTG",
        "region_modelMTG:phenotype_value"
      )
    )
    if (
      length(phenotype_index) != 1L ||
      length(mec_interaction_index) != 1L ||
      length(mtg_interaction_index) != 1L
    ) {
      stop("Repeated-model pathology contrasts are not identifiable")
    }
    contrast_matrix <- matrix(
      0,
      nrow = ncol(design),
      ncol = 4L,
      dimnames = list(
        coefficient_names,
        c("average_pathology", "DFC", "MEC", "MTG")
      )
    )
    contrast_matrix[phenotype_index, ] <- 1
    contrast_matrix[mec_interaction_index, "average_pathology"] <- 1 / 3
    contrast_matrix[mtg_interaction_index, "average_pathology"] <- 1 / 3
    contrast_matrix[mec_interaction_index, "MEC"] <- 1
    contrast_matrix[mtg_interaction_index, "MTG"] <- 1
    main_fit <- eBayes(
      contrasts.fit(base_fit, contrast_matrix),
      robust = TRUE
    )
    coefficient <- which(
      colnames(main_fit$coefficients) == "average_pathology"
    )
    main_effect_label <- "equal_region_average_pathology"

    contrast_se <- function(column) {
      statistic <- as.numeric(main_fit$t[, column])
      coefficient_value <- as.numeric(main_fit$coefficients[, column])
      ifelse(
        is.finite(statistic) & statistic != 0,
        abs(coefficient_value / statistic),
        NA_real_
      )
    }
    region_contrast_table <- data.frame(
      gene_metadata[keep, , drop = FALSE],
      beta_D_DFC = as.numeric(main_fit$coefficients[, "DFC"]),
      SE_D_DFC = contrast_se("DFC"),
      P_DFC = as.numeric(main_fit$p.value[, "DFC"]),
      FDR_DFC = p.adjust(
        as.numeric(main_fit$p.value[, "DFC"]),
        method = "BH"
      ),
      beta_D_MEC = as.numeric(main_fit$coefficients[, "MEC"]),
      SE_D_MEC = contrast_se("MEC"),
      P_MEC = as.numeric(main_fit$p.value[, "MEC"]),
      FDR_MEC = p.adjust(
        as.numeric(main_fit$p.value[, "MEC"]),
        method = "BH"
      ),
      beta_D_MTG = as.numeric(main_fit$coefficients[, "MTG"]),
      SE_D_MTG = contrast_se("MTG"),
      P_MTG = as.numeric(main_fit$p.value[, "MTG"]),
      FDR_MTG = p.adjust(
        as.numeric(main_fit$p.value[, "MTG"]),
        method = "BH"
      ),
      stringsAsFactors = FALSE
    )

    interaction_matrix <- matrix(
      0,
      nrow = ncol(design),
      ncol = 2L,
      dimnames = list(
        coefficient_names,
        c("MEC_minus_DFC", "MTG_minus_DFC")
      )
    )
    interaction_matrix[
      mec_interaction_index, "MEC_minus_DFC"
    ] <- 1
    interaction_matrix[
      mtg_interaction_index, "MTG_minus_DFC"
    ] <- 1
    interaction_fit <- eBayes(
      contrasts.fit(base_fit, interaction_matrix),
      robust = TRUE
    )
    interaction_top <- topTable(
      interaction_fit,
      coef = 1:2,
      number = Inf,
      sort.by = "none",
      adjust.method = "BH"
    )
    interaction_table <- cbind(
      gene_metadata[keep, , drop = FALSE],
      data.frame(
        F = interaction_top$F,
        P_joint_region_interaction = interaction_top$P.Value,
        FDR_joint_region_interaction = interaction_top$adj.P.Val,
        stringsAsFactors = FALSE
      )
    )
  } else {
    main_fit <- eBayes(base_fit, robust = TRUE)
    coefficient <- which(
      colnames(main_fit$coefficients) == "phenotype_value"
    )
  }
  if (length(coefficient) != 1L) {
    stop("Primary pathology contrast was not uniquely identified")
  }
  table <- topTable(
    main_fit,
    coef = coefficient,
    number = Inf,
    sort.by = "none",
    adjust.method = "BH"
  )
  moderated_t <- as.numeric(main_fit$t[, coefficient])
  moderated_se <- ifelse(
    is.finite(moderated_t) & moderated_t != 0,
    abs(
      as.numeric(main_fit$coefficients[, coefficient]) /
        moderated_t
    ),
    NA_real_
  )
  result <- cbind(
    gene_metadata[keep, , drop = FALSE],
    data.frame(
      class = settings$cell_class,
      phenotype = phenotype_name,
      analysis = settings$analysis,
      region = if (settings$analysis == "region") {
        settings$region
      } else {
        "DFC_MEC_MTG"
      },
      composition_adjusted = settings$composition_pcs,
      main_effect = main_effect_label,
      beta_D = as.numeric(main_fit$coefficients[, coefficient]),
      SE_D_moderated = moderated_se,
      t = moderated_t,
      AveExpr = table$AveExpr,
      P.Value = table$P.Value,
      FDR_within_class = table$adj.P.Val,
      B = table$B,
      n_samples = nrow(metadata),
      n_donors = length(unique(metadata$donor_id)),
      consensus_correlation = consensus_correlation,
      stringsAsFactors = FALSE
    )
  )

  result_path <- file.path(
    output_dir,
    sprintf(
      "%s__%s__gene_results.csv",
      settings$cell_class,
      phenotype_name
    )
  )
  sample_path <- file.path(
    output_dir,
    sprintf(
      "%s__%s__samples.csv",
      settings$cell_class,
      phenotype_name
    )
  )
  write_csv_atomic(result, result_path)
  region_contrast_path <- ""
  interaction_path <- ""
  if (!is.null(region_contrast_table)) {
    region_contrast_path <- file.path(
      output_dir,
      sprintf(
        "%s__%s__region_contrasts.csv",
        settings$cell_class,
        phenotype_name
      )
    )
    write_csv_atomic(region_contrast_table, region_contrast_path)
  }
  if (!is.null(interaction_table)) {
    interaction_path <- file.path(
      output_dir,
      sprintf(
        "%s__%s__joint_region_interaction.csv",
        settings$cell_class,
        phenotype_name
      )
    )
    write_csv_atomic(interaction_table, interaction_path)
  }
  write_csv_atomic(
    metadata[
      ,
      unique(c(
        "sample_id",
        "donor_id",
        "region",
        "phenotype_raw",
        "n_cells_numeric",
        design_record$composition_columns
      )),
      drop = FALSE
    ],
    sample_path
  )
  write_csv_atomic(
    cbind(
      data.frame(sample_id = metadata$sample_id),
      as.data.frame(design)
    ),
    file.path(
      output_dir,
      sprintf(
        "%s__%s__design_matrix.csv",
        settings$cell_class,
        phenotype_name
      )
    )
  )

  data.frame(
    class = settings$cell_class,
    phenotype = phenotype_name,
    status = "complete",
    analysis = settings$analysis,
    region = if (settings$analysis == "region") {
      settings$region
    } else {
      "DFC_MEC_MTG"
    },
    composition_adjusted = settings$composition_pcs,
    composition_columns = paste(
      design_record$composition_columns,
      collapse = ";"
    ),
    model_formula = design_record$formula,
    n_samples = nrow(metadata),
    n_donors = length(unique(metadata$donor_id)),
    phenotype_min = min(metadata$phenotype_raw),
    phenotype_max = max(metadata$phenotype_raw),
    genes_tested = nrow(result),
    design_rank = qr(design)$rank,
    design_columns = ncol(design),
    design_condition_number = kappa(design, exact = TRUE),
    minimum_gene_samples = minimum_samples,
    consensus_correlation = consensus_correlation,
    main_effect = main_effect_label,
    fdr05_genes = sum(
      result$FDR_within_class < 0.05,
      na.rm = TRUE
    ),
    joint_interaction_fdr05_genes = if (
      is.null(interaction_table)
    ) {
      NA_integer_
    } else {
      sum(
        interaction_table$FDR_joint_region_interaction < 0.05,
        na.rm = TRUE
      )
    },
    result_file = basename(result_path),
    region_contrast_file = if (nzchar(region_contrast_path)) {
      basename(region_contrast_path)
    } else {
      ""
    },
    joint_interaction_file = if (nzchar(interaction_path)) {
      basename(interaction_path)
    } else {
      ""
    },
    sample_file = basename(sample_path),
    stringsAsFactors = FALSE
  )
}

settings <- parse_args(commandArgs(trailingOnly = TRUE))
input <- normalizePath(settings$input, mustWork = TRUE)
output_dir <- normalizePath(settings$output_dir, mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading ", input)
sce <- readH5AD(input, use_hdf5 = FALSE, reader = "R")
if (!"X" %in% assayNames(sce)) {
  stop("X assay missing")
}
counts <- assay(sce, "X")
if (!inherits(counts, "sparseMatrix")) {
  counts <- as(counts, "dgCMatrix")
}
if (nrow(counts) != 36601L) {
  stop("Unexpected gene count")
}
if (
  any(counts@x < 0) ||
  any(abs(counts@x - round(counts@x)) > 1e-6)
) {
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
if (anyDuplicated(gene_metadata$gene_id)) {
  stop("Duplicated gene IDs")
}

manifest_rows <- list()
for (phenotype_name in settings$phenotypes) {
  result_file <- file.path(
    output_dir,
    sprintf(
      "%s__%s__gene_results.csv",
      settings$cell_class,
      phenotype_name
    )
  )
  if (file.exists(result_file) && !settings$overwrite) {
    message("Skipping existing ", basename(result_file))
    next
  }
  metadata <- prepare_metadata(
    colData(sce),
    phenotype_name,
    settings
  )
  message(
    "Fitting ",
    settings$cell_class,
    " / ",
    phenotype_name,
    " / ",
    settings$analysis,
    if (settings$analysis == "region") {
      paste0(" / ", settings$region)
    } else {
      ""
    },
    if (settings$composition_pcs) " / composition PCs" else ""
  )
  manifest_rows[[phenotype_name]] <- fit_one(
    counts,
    gene_metadata,
    metadata,
    phenotype_name,
    settings,
    output_dir
  )
}

if (length(manifest_rows)) {
  manifest <- do.call(rbind, manifest_rows)
  write_csv_atomic(
    manifest,
    file.path(output_dir, "model_manifest.csv")
  )
}
sink(file.path(output_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()
message("SEA-AD donor-region model pipeline complete: ", output_dir)
