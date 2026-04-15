# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains the source code for **sigQC**, an R package for quality control of gene signatures across genomic datasets. The methodology is published in Nature Protocols (Dhawan et al., 2019). The project directory name "pysigQC" reflects an intent to create a Python port of this R package.

The R source lives under `sigQC-master/`.

## R Package Commands

```bash
# Install from source
R CMD INSTALL sigQC-master/

# Or via devtools in R
Rscript -e 'devtools::install("sigQC-master")'

# Check package (build + tests + examples)
R CMD check sigQC-master/

# Run tests (testthat v3)
Rscript -e 'devtools::test("sigQC-master")'

# Build documentation (roxygen2)
Rscript -e 'devtools::document("sigQC-master")'
```

## Architecture

### Single entry point

The package exports exactly one function: `make_all_plots()` (in `R/make_all_plots.R`). It orchestrates the entire QC pipeline by sequentially calling internal `*_loc` analysis functions, accumulating radar plot summary values as it goes.

### Analysis modules (all internal, in `R/`)

Each module follows the same interface pattern: `eval_*_loc(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets, ...)` and produces PDF plots plus text/CSV output tables.

| Module | Purpose |
|---|---|
| `eval_expr_loc.R` | Expression levels: proportion expressed, NA analysis, density plots |
| `eval_var_loc.R` | Variability: mean vs SD, coefficient of variation, skewness |
| `eval_compactness_loc.R` | Compactness: autocorrelation heatmaps, rank products across datasets |
| `eval_struct_loc.R` | Structure: hierarchical clustering heatmaps, biclustering on binarized data |
| `compare_metrics_loc.R` | Scoring comparison: mean, median, PCA1 correlation, scree plots |
| `eval_stan_loc.R` | Standardization effects: raw vs z-transformed score correlations |
| `make_radar_chart_loc.R` | Summary radar chart integrating all computed metrics |
| `sigsQcNegativeControl.R` | Negative/permutation controls via bootstrap resampling |

### Supporting utilities

- `compute_without_plots.R` — runs computations without generating visualizations
- `draw.heatmaps.R` — heatmap rendering helpers
- `boxplot.matrix2.R` — extended boxplot for matrices
- `grab_grob.R` — grid graphics object capture

### Data flow

1. `make_all_plots()` validates inputs (signatures list, expression matrices, names)
2. Each `eval_*_loc` function runs independently, writing PDFs and tables to `out_dir`
3. Each module returns metric values that accumulate into `radar_plot_values`
4. `make_radar_chart_loc` produces the final summary radar chart from accumulated metrics
5. Optionally, `sigsQcNegativeControl` runs bootstrap resampling for null distributions

### Input requirements

- Gene signatures: list of character vectors (gene IDs must match expression matrix rownames)
- Expression matrices: list of numeric matrices (genes × samples), should be normalized, batch-corrected, and log-transformed prior to use
- Minimum 2 genes per signature, 2 samples per dataset

### Key dependencies

Bioconductor: `ComplexHeatmap`, `GSVA`. CRAN: `biclust`, `fmsb`, `moments`, `mclust`, `circlize`, `gplots`. Optional: `RankProd` for rank product analysis across datasets.

