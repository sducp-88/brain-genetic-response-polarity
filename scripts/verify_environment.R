required_packages <- c(
  "dreamlet",
  "variancePartition",
  "SingleCellExperiment",
  "edgeR",
  "limma",
  "DelayedArray",
  "HDF5Array",
  "zellkonverter",
  "coloc",
  "susieR",
  "data.table",
  "renv"
)

package_status <- data.frame(
  package = required_packages,
  installed = vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  ),
  version = vapply(
    required_packages,
    function(package_name) {
      if (requireNamespace(package_name, quietly = TRUE)) {
        as.character(packageVersion(package_name))
      } else {
        NA_character_
      }
    },
    character(1)
  ),
  row.names = NULL
)

print(package_status, row.names = FALSE)

if (!all(package_status$installed)) {
  missing_packages <- package_status$package[!package_status$installed]
  stop(
    "Missing R packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

cat("\nR environment verification passed.\n")
sessionInfo()
