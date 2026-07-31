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
    result_dir = NULL,
    output_dir = NULL,
    min_cells = 20L,
    permutations = 20L,
    workers = 2L,
    seed = 20260729L
  )
  i <- 1L
  while (i <= length(x)) {
    key <- gsub("-", "_", sub("^--", "", x[[i]]))
    if (i == length(x)) stop("Missing value for --", key)
    out[[key]] <- x[[i + 1L]]
    i <- i + 2L
  }
  for (key in c("input_dir", "result_dir", "output_dir")) {
    if (is.null(out[[key]])) stop("--", gsub("_", "-", key), " is required")
  }
  out$min_cells <- as.integer(out$min_cells)
  out$permutations <- as.integer(out$permutations)
  out$workers <- as.integer(out$workers)
  out$seed <- as.integer(out$seed)
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
    stop("Only continuous pathology phenotypes are supported")
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

lambda_gc <- function(p) {
  p <- p[is.finite(p) & p >= 0 & p <= 1]
  p <- pmax(p, .Machine$double.xmin)
  median(qchisq(p, df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1)
}

settings <- parse_args(commandArgs(trailingOnly = TRUE))
input_dir <- normalizePath(settings$input_dir, mustWork = TRUE)
result_dir <- normalizePath(settings$result_dir, mustWork = TRUE)
output_dir <- normalizePath(settings$output_dir, mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(settings$seed)

classes <- c("Astro", "EN", "Endo", "IN", "Immune", "Mural", "OPC", "Oligo")
phenotypes <- c("CERAD", "Braak")
all_permutations <- list()
summary_rows <- list()

for (class_index in seq_along(classes)) {
  class_name <- classes[[class_index]]
  class_complete <- all(vapply(phenotypes, function(phenotype_name) {
    checkpoint_path <- file.path(
      output_dir,
      sprintf("%s__%s__permutations.csv", class_name, phenotype_name)
    )
    summary_checkpoint_path <- file.path(
      output_dir,
      sprintf("%s__%s__summary.csv", class_name, phenotype_name)
    )
    if (!file.exists(checkpoint_path) || !file.exists(summary_checkpoint_path)) {
      return(FALSE)
    }
    nrow(read.csv(checkpoint_path, stringsAsFactors = FALSE)) ==
      settings$permutations
  }, logical(1)))
  if (class_complete) {
    for (phenotype_name in phenotypes) {
      checkpoint_path <- file.path(
        output_dir,
        sprintf("%s__%s__permutations.csv", class_name, phenotype_name)
      )
      summary_checkpoint_path <- file.path(
        output_dir,
        sprintf("%s__%s__summary.csv", class_name, phenotype_name)
      )
      all_permutations[[length(all_permutations) + 1L]] <- read.csv(
        checkpoint_path,
        stringsAsFactors = FALSE
      )
      summary_rows[[length(summary_rows) + 1L]] <- read.csv(
        summary_checkpoint_path,
        stringsAsFactors = FALSE
      )
    }
    message("Skipping completed class ", class_name)
    next
  }
  input_file <- file.path(
    input_dir, sprintf("RADC_%s_pseudobulk_counts.h5ad", class_name)
  )
  message("Reading ", input_file)
  sce <- readH5AD(input_file, use_hdf5 = FALSE, reader = "R")
  counts <- assay(sce, "X")
  if (!inherits(counts, "sparseMatrix")) counts <- as(counts, "dgCMatrix")
  storage.mode(counts@x) <- "double"

  for (phenotype_index in seq_along(phenotypes)) {
    phenotype_name <- phenotypes[[phenotype_index]]
    checkpoint_path <- file.path(
      output_dir,
      sprintf("%s__%s__permutations.csv", class_name, phenotype_name)
    )
    summary_checkpoint_path <- file.path(
      output_dir,
      sprintf("%s__%s__summary.csv", class_name, phenotype_name)
    )
    if (file.exists(checkpoint_path) && file.exists(summary_checkpoint_path)) {
      permutation_frame <- read.csv(checkpoint_path, stringsAsFactors = FALSE)
      summary_row <- read.csv(summary_checkpoint_path, stringsAsFactors = FALSE)
      if (
        nrow(permutation_frame) == settings$permutations &&
        nrow(summary_row) == 1L
      ) {
        message("Skipping completed ", class_name, " / ", phenotype_name)
        all_permutations[[length(all_permutations) + 1L]] <- permutation_frame
        summary_rows[[length(summary_rows) + 1L]] <- summary_row
        next
      }
    }
    message("Calibrating ", class_name, " / ", phenotype_name)
    metadata <- prepare_metadata(colData(sce), phenotype_name, settings$min_cells)
    columns <- match(metadata$sample_id, colnames(counts))
    model_counts <- counts[, columns, drop = FALSE]
    design <- model.matrix(
      ~ phenotype_value + age_z + sex_model + pmi_z + ancestry_model +
        log10_n_cells_z + nonAD_dx,
      data = metadata
    )
    covariate_design <- model.matrix(
      ~ age_z + sex_model + pmi_z + ancestry_model + log10_n_cells_z +
        nonAD_dx,
      data = metadata
    )
    dge <- DGEList(counts = model_counts)
    minimum_samples <- max(10L, ceiling(0.10 * ncol(dge)))
    keep_gene <- rowSums(cpm(dge) >= 1) >= minimum_samples
    dge <- dge[keep_gene, , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge, method = "TMM")
    voom_object <- voom(dge, design = design, plot = FALSE, normalize.method = "none")

    phenotype_fit <- lm.fit(covariate_design, metadata$phenotype_value)
    fitted_value <- as.numeric(phenotype_fit$fitted.values)
    residual_value <- as.numeric(phenotype_fit$residuals)
    model_seed <- settings$seed + class_index * 1000L + phenotype_index * 100L
    set.seed(model_seed)
    residual_orders <- lapply(
      seq_len(settings$permutations),
      function(i) sample.int(length(residual_value), replace = FALSE)
    )
    permutation_results <- parallel::mclapply(
      seq_len(settings$permutations),
      function(permutation_index) {
        perm_metadata <- metadata
        perm_metadata$phenotype_value <- fitted_value +
          residual_value[residual_orders[[permutation_index]]]
        perm_design <- model.matrix(
          ~ phenotype_value + age_z + sex_model + pmi_z + ancestry_model +
            log10_n_cells_z + nonAD_dx,
          data = perm_metadata
        )
        fit <- eBayes(lmFit(voom_object, perm_design), robust = TRUE)
        coefficient <- which(colnames(perm_design) == "phenotype_value")
        p <- fit$p.value[, coefficient]
        data.frame(
          class = class_name,
          phenotype = phenotype_name,
          permutation = permutation_index,
          lambda_gc = lambda_gc(p),
          p_lt_0_05 = sum(p < 0.05),
          p_lt_0_001 = sum(p < 0.001),
          fdr_lt_0_05 = sum(p.adjust(p, method = "BH") < 0.05),
          minimum_p = min(p),
          stringsAsFactors = FALSE
        )
      },
      mc.cores = min(settings$workers, settings$permutations),
      mc.preschedule = FALSE,
      mc.set.seed = FALSE
    )
    permutation_frame <- do.call(rbind, permutation_results)

    actual_path <- file.path(
      result_dir,
      sprintf("%s__%s__gene_results.csv", class_name, phenotype_name)
    )
    actual <- read.csv(actual_path, stringsAsFactors = FALSE)
    actual_lambda <- lambda_gc(actual$P.Value)
    summary_row <- data.frame(
      class = class_name,
      phenotype = phenotype_name,
      n_donors = nrow(metadata),
      genes = nrow(actual),
      actual_lambda_gc = actual_lambda,
      permutation_lambda_median = median(permutation_frame$lambda_gc),
      permutation_lambda_q025 = quantile(permutation_frame$lambda_gc, 0.025),
      permutation_lambda_q975 = quantile(permutation_frame$lambda_gc, 0.975),
      empirical_p_lambda = (
        1 + sum(permutation_frame$lambda_gc >= actual_lambda)
      ) / (1 + settings$permutations),
      actual_fdr05 = sum(actual$FDR_within_class < 0.05),
      permutation_fdr05_median = median(permutation_frame$fdr_lt_0_05),
      permutation_fdr05_max = max(permutation_frame$fdr_lt_0_05),
      permutations = settings$permutations,
      stringsAsFactors = FALSE
    )
    write.csv(permutation_frame, checkpoint_path, row.names = FALSE)
    write.csv(summary_row, summary_checkpoint_path, row.names = FALSE)
    all_permutations[[length(all_permutations) + 1L]] <- permutation_frame
    summary_rows[[length(summary_rows) + 1L]] <- summary_row
    rm(model_counts, dge, voom_object)
    gc()
  }
  rm(sce, counts)
  gc()
}

permutations <- do.call(rbind, all_permutations)
summary <- do.call(rbind, summary_rows)
summary$empirical_p_fdr_count <- vapply(seq_len(nrow(summary)), function(i) {
  selected <- permutations$class == summary$class[[i]] &
    permutations$phenotype == summary$phenotype[[i]]
  (
    1 + sum(
      permutations$fdr_lt_0_05[selected] >= summary$actual_fdr05[[i]]
    )
  ) / (1 + sum(selected))
}, numeric(1))
summary$actual_lambda_above_all_permutations <- vapply(
  seq_len(nrow(summary)),
  function(i) {
    selected <- permutations$class == summary$class[[i]] &
      permutations$phenotype == summary$phenotype[[i]]
    summary$actual_lambda_gc[[i]] > max(permutations$lambda_gc[selected])
  },
  logical(1)
)
summary$actual_fdr_count_above_all_permutations <- vapply(
  seq_len(nrow(summary)),
  function(i) {
    selected <- permutations$class == summary$class[[i]] &
      permutations$phenotype == summary$phenotype[[i]]
    summary$actual_fdr05[[i]] > max(permutations$fdr_lt_0_05[selected])
  },
  logical(1)
)
write.csv(
  permutations,
  file.path(output_dir, "permutation_diagnostics.csv"),
  row.names = FALSE
)
write.csv(
  summary,
  file.path(output_dir, "permutation_calibration_summary.csv"),
  row.names = FALSE
)

png(
  file.path(output_dir, "lambda_permutation_calibration.png"),
  width = 2100,
  height = 1300,
  res = 180
)
old_par <- par(mar = c(10, 5, 3, 1))
model_labels <- paste(summary$class, summary$phenotype, sep = "\n")
permutations$model <- factor(
  paste(permutations$class, permutations$phenotype, sep = "\n"),
  levels = model_labels
)
ylim <- range(
  c(
    permutations$lambda_gc,
    summary$actual_lambda_gc
  ),
  finite = TRUE
)
boxplot(
  lambda_gc ~ model,
  data = permutations,
  las = 2,
  ylab = "Lambda GC",
  xlab = "",
  main = "Freedman-Lane phenotype-residual permutation calibration",
  ylim = ylim,
  outline = FALSE
)
points(
  seq_len(nrow(summary)),
  summary$actual_lambda_gc,
  pch = 19,
  col = "red"
)
abline(h = 1, lty = 2, col = "grey40")
legend("topleft", legend = "Observed", pch = 19, col = "red", bty = "n")
par(old_par)
dev.off()

report <- c(
  "# RADC continuous-pathology permutation calibration",
  "",
  paste0("- Models: ", nrow(summary), "/16."),
  paste0("- Permutations per model: ", settings$permutations, "."),
  paste0(
    "- Models with observed lambda above every permutation: ",
    sum(summary$actual_lambda_above_all_permutations),
    "/16."
  ),
  paste0(
    "- Models with observed FDR-hit count above every permutation: ",
    sum(summary$actual_fdr_count_above_all_permutations),
    "/16."
  ),
  "",
  "Freedman-Lane residual permutation preserves phenotype-covariate structure.",
  "This is an initial calibration audit; publication-grade tail calibration",
  "requires more permutations for prioritized models."
)
writeLines(report, file.path(output_dir, "PERMUTATION_CALIBRATION_REPORT.md"))
sink(file.path(output_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()
message("Permutation calibration complete: ", output_dir)
