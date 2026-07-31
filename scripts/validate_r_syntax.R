#!/usr/bin/env Rscript

files <- list.files(
  "scripts",
  pattern = "[.]R$",
  recursive = TRUE,
  full.names = TRUE
)

invisible(lapply(files, parse))
cat("R files parsed:", length(files), "\n")
