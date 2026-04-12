# testthat.R — Test runner
# Run with: Rscript tests/testthat.R  (from sigQC-master/)
# Or: Rscript -e 'source("tests/testthat.R")'

library(testthat)

# Load required dependencies
library(mclust)
library(GSVA)
library(biclust)

# Source all refactored compute_*() functions
script_args <- commandArgs(trailingOnly = FALSE)
file_flag   <- grep("^--file=", script_args, value = TRUE)
if (length(file_flag) > 0) {
  script_path  <- normalizePath(sub("^--file=", "", file_flag[1]))
  refactored_dir <- file.path(dirname(dirname(script_path)), "R_refactored")
} else {
  refactored_dir <- file.path(dirname(dirname(normalizePath("."))), "R_refactored")
}
if (!dir.exists(refactored_dir)) {
  refactored_dir <- "R_refactored"
}

for (f in list.files(refactored_dir, pattern = "[.]R$", full.names = TRUE)) {
  source(f)
}

test_dir("tests/testthat", reporter = "summary")
