# testthat.R — Test runner
# Run with: Rscript tests/testthat.R  (from sigQC-master/)
# Or: Rscript -e 'source("tests/testthat.R")'

library(testthat)

# Load required dependencies
library(mclust)
library(GSVA)
library(biclust)

# Source all refactored compute_*() functions
refactored_dir <- file.path(dirname(dirname(sys.frame(1)$ofile %||% ".")), "R_refactored")
if (!dir.exists(refactored_dir)) {
  refactored_dir <- "R_refactored"
}

for (f in list.files(refactored_dir, pattern = "[.]R$", full.names = TRUE)) {
  source(f)
}

test_dir("tests/testthat", reporter = "summary")
