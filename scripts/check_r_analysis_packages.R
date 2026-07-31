packages <- c(
  "dreamlet",
  "variancePartition",
  "edgeR",
  "limma",
  "SingleCellExperiment",
  "zellkonverter",
  "HDF5Array",
  "Matrix"
)

for (package in packages) {
  available <- requireNamespace(package, quietly = TRUE)
  version <- if (available) as.character(utils::packageVersion(package)) else ""
  cat(package, available, version, "\n")
}
