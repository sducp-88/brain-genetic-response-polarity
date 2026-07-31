suppressPackageStartupMessages({
  library(dreamlet)
  library(SingleCellExperiment)
})

set.seed(20260725)

n_genes <- 300L
n_donors <- 16L
cells_per_type <- 50L
cell_types <- c("Neuron", "Microglia")

donor_ids <- sprintf("D%02d", seq_len(n_donors))
donor_diagnosis <- factor(
  rep(c("Control", "Case"), each = n_donors / 2L),
  levels = c("Control", "Case")
)

cell_metadata <- do.call(
  rbind,
  lapply(
    seq_along(donor_ids),
    function(donor_index) {
      data.frame(
        donor = donor_ids[[donor_index]],
        diagnosis = donor_diagnosis[[donor_index]],
        cell_type = rep(cell_types, each = cells_per_type),
        stringsAsFactors = FALSE
      )
    }
  )
)

cell_metadata$diagnosis <- factor(
  cell_metadata$diagnosis,
  levels = c("Control", "Case")
)

n_cells <- nrow(cell_metadata)
counts_matrix <- matrix(
  rpois(n_genes * n_cells, lambda = 1.5),
  nrow = n_genes,
  ncol = n_cells,
  dimnames = list(
    sprintf("GENE%03d", seq_len(n_genes)),
    sprintf("CELL%05d", seq_len(n_cells))
  )
)

signal_genes <- sprintf("GENE%03d", seq_len(12L))
signal_cells <- which(
  cell_metadata$diagnosis == "Case" &
    cell_metadata$cell_type == "Neuron"
)

counts_matrix[signal_genes, signal_cells] <-
  counts_matrix[signal_genes, signal_cells] +
  matrix(
    rpois(length(signal_genes) * length(signal_cells), lambda = 2.5),
    nrow = length(signal_genes)
  )

single_cell_object <- SingleCellExperiment(
  assays = list(counts = counts_matrix),
  colData = S4Vectors::DataFrame(cell_metadata)
)

pseudobulk_object <- aggregateToPseudoBulk(
  single_cell_object,
  assay = "counts",
  cluster_id = "cell_type",
  sample_id = "donor",
  verbose = FALSE
)

processed_object <- processAssays(
  pseudobulk_object,
  ~ diagnosis,
  min.count = 1
)

model_result <- dreamlet(
  processed_object,
  ~ diagnosis
)

coefficient_name <- "diagnosisCase"

result_table <- as.data.frame(
  topTable(
    model_result,
    coef = coefficient_name,
    number = Inf
  )
)

neuron_results <- result_table[
  result_table$assay == "Neuron" &
    result_table$ID %in% signal_genes,
  ,
  drop = FALSE
]

if (nrow(neuron_results) != length(signal_genes)) {
  stop(
    "Not all simulated signal genes were retained in the neuron model.",
    call. = FALSE
  )
}

if (mean(neuron_results$logFC > 0) < 0.9) {
  stop(
    "Dreamlet did not recover the expected positive direction.",
    call. = FALSE
  )
}

script_argument <- commandArgs(trailingOnly = TRUE)
project_directory <- if (length(script_argument) >= 1L) {
  normalizePath(script_argument[[1L]], mustWork = TRUE)
} else {
  normalizePath(file.path(getwd()), mustWork = TRUE)
}

output_directory <- file.path(project_directory, "outputs")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

output_file <- file.path(
  output_directory,
  "smoke_test_dreamlet_results.csv"
)

utils::write.csv(
  result_table,
  output_file,
  row.names = FALSE
)

cat("Dreamlet smoke test passed.\n")
cat("Donors:", n_donors, "\n")
cat("Cells:", n_cells, "\n")
cat("Genes:", n_genes, "\n")
cat(
  "Recovered positive signal genes:",
  sum(neuron_results$logFC > 0),
  "of",
  length(signal_genes),
  "\n"
)
cat("Output:", output_file, "\n")
