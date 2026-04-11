# sigQC Project Guide

Comprehensive documentation for the sigQC audit, refactoring, testing, and Python port.

---

## Table of Contents

1. [Project Overview](#project-overview)
2. [Directory Structure](#directory-structure)
3. [Code Versions](#code-versions)
4. [Bugs Found (AUDIT.md)](#bugs-found)
5. [Test Fixtures](#test-fixtures)
6. [Running R Tests](#running-r-tests)
7. [Running Python Tests](#running-python-tests)
8. [Cross-Validation: R vs Python](#cross-validation-r-vs-python)
9. [SLURM Cluster Execution](#slurm-cluster-execution)
10. [Analyzing Discrepancies](#analyzing-discrepancies)
11. [Python Package Usage](#python-package-usage)
12. [Key R/Python Behavioral Differences](#key-rpython-behavioral-differences)

---

## Project Overview

**sigQC** is an R package for quality control of gene signatures across genomic datasets (Dhawan et al., Nature Protocols 2019). This project:

1. **Audited** the R source — found 8 bugs and 7 theoretical issues (see `sigQC-master/AUDIT.md`)
2. **Corrected** the bugs in `R_corrected/`
3. **Refactored** the code into testable compute/plot separation in `R_refactored/`
4. **Optimized** the negative control module in `R_optimized/`
5. **Ported** all compute functions to Python (`pysigqc/`) with corrections and optimizations
6. **Cross-validated** Python outputs against R reference outputs on shared fixtures

---

## Directory Structure

```
pysigQC/
├── CLAUDE.md                  # Claude Code project instructions
├── GUIDE.md                   # This file
├── pyproject.toml             # Python package configuration
│
├── pysigqc/                   # Python port (corrected + optimized)
│   ├── __init__.py
│   ├── utils.py               # Shared helpers (gene_intersection, z_transform)
│   ├── eval_var.py            # Variability metrics (6 radar values)
│   ├── eval_expr.py           # Expression metrics (2 radar values)
│   ├── eval_compactness.py    # Compactness metrics (1 radar value)
│   ├── eval_stan.py           # Standardization metrics (1 radar value)
│   ├── compare_metrics.py     # Scoring comparison metrics (4 radar values)
│   ├── eval_struct.py         # Structure analysis (clustering, biclustering)
│   ├── radar_chart.py         # Radar chart assembly (14 metrics)
│   └── negative_control.py    # Negative/permutation controls (optimized)
│
├── tests/                     # Python tests
│   ├── conftest.py            # Shared pytest fixtures
│   ├── test_eval_var.py
│   ├── test_eval_expr.py
│   ├── test_eval_compactness.py
│   ├── test_eval_stan.py
│   ├── test_compare_metrics.py
│   ├── test_eval_struct.py
│   ├── test_radar_chart.py
│   ├── test_negative_control.py
│   └── test_medium_cross_validation.py  # Medium fixture cross-validation
│
└── sigQC-master/              # R package and all R code versions
    ├── DESCRIPTION / NAMESPACE         # Point to R/ (original package installable)
    ├── AUDIT.md                        # Bug and issue documentation
    ├── R/                              # ORIGINAL (untouched, bugs included)
    ├── R_corrected/                    # Bug fixes only (8 bugs fixed)
    ├── R_refactored/                   # Corrected + compute/plot split + testable
    ├── R_optimized/                    # Optimized sigsQcNegativeControl.R only
    ├── tests/
    │   ├── testthat.R                  # R test runner
    │   ├── fixtures/
    │   │   ├── fixture_generator.R             # Small fixture (set.seed(42))
    │   │   ├── fixture_medium_generator.R      # Medium fixture (set.seed(2024))
    │   │   ├── generate_reference_outputs.R    # Small R reference output generator
    │   │   ├── generate_medium_reference.R     # Medium R reference output generator
    │   │   ├── fixture_small.rds / fixture_*.csv
    │   │   ├── fixture_medium.rds / fixture_medium_*.csv
    │   │   ├── reference_outputs/              # 43 small reference CSVs
    │   │   └── reference_outputs_medium/       # 76 medium reference CSVs
    │   └── testthat/
    │       ├── test-compute_var.R
    │       ├── test-compute_expr.R
    │       ├── test-compute_compactness.R
    │       ├── test-compute_stan.R
    │       ├── test-compute_metrics.R
    │       ├── test-compute_struct.R
    │       ├── test-compute_radar.R
    │       ├── test-negative_control.R
    │       └── test-integration.R
    ├── slurm/
    │   ├── run_r_tests.slurm               # SLURM: R fixture generation + testthat
    │   ├── run_py_tests.slurm              # SLURM: Python pytest suite
    │   ├── run_brca_cross_validation.slurm # SLURM: R + Python on BRCA data
    │   ├── run_original.slurm              # SLURM: original R pipeline
    │   └── run_corrected.slurm             # SLURM: corrected R pipeline
    └── input_files/                        # BRCA paper data (not in git)
```

---

## Code Versions

| Version | Location | Purpose | Bugs Fixed? | Testable? |
|---------|----------|---------|-------------|-----------|
| **Original** | `R/` | Untouched source as published | No | No (compute + plot interleaved) |
| **Corrected** | `R_corrected/` | Bug fixes only, minimal changes | Yes (8 bugs) | No (same structure as original) |
| **Refactored** | `R_refactored/` | Corrected + `compute_*()`/`plot_*()` split | Yes | Yes (compute functions return data) |
| **Optimized** | `R_optimized/` | Optimized negative control only | Yes | Yes |
| **Python** | `pysigqc/` | Full port, corrected + optimized | Yes | Yes |

The `DESCRIPTION` and `NAMESPACE` files still point to `R/`, so `R CMD INSTALL sigQC-master/` installs the original package.

---

## Bugs Found

See `sigQC-master/AUDIT.md` for full details. Summary of the 8 confirmed bugs:

| ID | Severity | File | Issue |
|----|----------|------|-------|
| BUG-1 | CRITICAL | `sigsQcNegativeControl.R` | `runif()` instead of `sample()` for gene sampling |
| BUG-2 | CRITICAL | `eval_stan_loc.R`, `eval_struct_loc.R` | Division by zero for zero-variance genes |
| BUG-3 | CRITICAL | `make_all_plots.R` | Index bug when removing short signatures |
| BUG-4 | MEDIUM | `compare_metrics_loc.R` | Mclust E/V model labels swapped |
| BUG-5 | MEDIUM | `eval_struct_loc.R` | `biclust()` called with wrong matrix dimensions |
| BUG-6 | LOW | `eval_expr_loc.R` | Threshold computed on all genes instead of signature genes |
| BUG-7 | LOW | `eval_var_loc.R` | CV computed as `mean/sd` instead of `sd/mean` |
| BUG-8 | LOW | `eval_struct_loc.R` | Missing genes not padded with NA for asymmetric datasets |

All bugs are fixed in `R_corrected/`, `R_refactored/`, and `pysigqc/`.

---

## Test Fixtures

### Small Fixture (primary)

- **Generator:** `sigQC-master/tests/fixtures/fixture_generator.R`
- **Seed:** `set.seed(42)`
- **Dimensions:** 10 genes, 10 samples, 2 datasets, 2 signatures
- **Edge cases:** NA values, one zero-variance gene (BUG-2)
- **Outputs:** `fixture_small.rds`, `fixture_dataset_A.csv`, `fixture_dataset_B.csv`, `fixture_signatures.csv`

### Medium Fixture (stress test)

- **Generator:** `sigQC-master/tests/fixtures/fixture_medium_generator.R`
- **Seed:** `set.seed(2024)`
- **Dimensions:** 100 genes, 50 samples, 3 datasets, 5 signatures
- **Edge cases:**
  - Scattered NAs (5 random positions in dataset A)
  - Block NAs (gene 15, last 10 samples in dataset B)
  - Full-row NA (gene 50 in dataset C)
  - Zero-variance gene (gene 20 in dataset B)
  - Near-zero-variance gene (gene 21 in dataset B, sd ~1e-8)
  - Missing genes (dataset C has only 90 of 100 genes)
  - Batch effects between datasets A, B, C
  - 5 signature types: tight co-expression, loose co-expression, mixed quality, high variance, random baseline

### Reference Outputs

Generated by running R `compute_*()` functions on the fixtures:

- `reference_outputs/` — 43 CSV files from small fixture
- `reference_outputs_medium/` — 76 CSV files from medium fixture

These are the ground truth for Python cross-validation.

---

## Running R Tests

### Prerequisites

R 4.3+ with packages: `testthat`, `mclust`, `GSVA`, `biclust`, `moments`, `fmsb`, `circlize`, `ComplexHeatmap`, `gplots`

### Local Execution

```bash
cd sigQC-master

# Step 1: Generate fixtures (deterministic)
Rscript tests/fixtures/fixture_generator.R
Rscript tests/fixtures/fixture_medium_generator.R

# Step 2: Generate reference outputs
Rscript tests/fixtures/generate_reference_outputs.R
Rscript tests/fixtures/generate_medium_reference.R

# Step 3: Run testthat
Rscript tests/testthat.R
```

The test runner sources all files from `R_refactored/` and runs 9 test files covering all compute modules + integration.

### On SLURM Cluster

```bash
cd sigQC-master/slurm
sbatch run_r_tests.slurm
```

This creates the conda environment `sigqc_env` (if needed), installs all R dependencies, generates fixtures, generates reference outputs, and runs testthat. Output goes to `slurm/out/`.

---

## Running Python Tests

### Prerequisites

```bash
# From project root (pysigQC/)
pip install -e ".[dev]"
```

This installs: `numpy`, `scipy`, `pandas`, `scikit-learn`, `matplotlib`, `seaborn`, `joblib`, `pytest`, `pytest-cov`.

### Local Execution

```bash
# From project root (pysigQC/)

# All tests (unit + cross-validation)
pytest tests/ -v

# Just unit tests (no R reference needed)
pytest tests/ -v -k "not cross_validation and not medium"

# Just cross-validation against R references
pytest tests/ -v -k "cross_validation or medium"

# Single module
pytest tests/test_eval_var.py -v
```

**Important:** The cross-validation tests (`test_*_cross_validation` functions and `test_medium_cross_validation.py`) require R reference outputs to exist. Run the R fixture/reference generation first (Step 2 above). If reference files are missing, these tests are skipped (not failed).

### On SLURM Cluster

```bash
cd sigQC-master/slurm
sbatch run_py_tests.slurm
```

Creates `sigqc_py_env` conda environment with Python 3.11, installs pysigqc in dev mode, and runs the full pytest suite. Checks that fixtures exist first.

---

## Cross-Validation: R vs Python

Cross-validation verifies that Python produces the same numerical results as the corrected R code on identical inputs.

### How It Works

1. R fixture generators create deterministic CSV datasets + signatures
2. R `generate_reference_outputs.R` runs all `compute_*()` functions and saves radar values as CSVs
3. Python tests load the same CSV inputs, run Python `compute_*()`, and compare against R CSVs

### Tolerance Levels

| Metric type | Tolerance | Reason |
|-------------|-----------|--------|
| Most metrics | `rtol=1e-4` | Floating-point differences between R/Python implementations |
| PCA-related (`rho_pca1_*`, `prop_pca1_var`) | `abs()` comparison, `rtol=0.05` | PCA sign ambiguity (eigenvector direction is arbitrary) |
| Mclust/GaussianMixture metrics | `rtol=0.05` | Different EM implementations |
| NaN values | Match if both NaN | NA propagation should be identical |

### Running BRCA Cross-Validation

For full validation on the original paper data:

```bash
cd sigQC-master/slurm
sbatch run_brca_cross_validation.slurm
```

This runs both R and Python on the BRCA data in `input_files/`, compares all 14 radar metrics, and reports mismatches to `brca_cross_validation/py_output/brca_mismatches.csv`.

**Prerequisite:** BRCA data must be in `sigQC-master/input_files/` with expression matrices as `exprs*.csv` and signatures as `sig*.csv`.

---

## Analyzing Discrepancies

When R and Python values differ, use this systematic approach:

### Step 1: Identify the Metric

Check which radar metric disagrees. The 14 metrics are:

| Module | Metrics |
|--------|---------|
| `eval_var` | `sd_median_ratio`, `abs_skewness_ratio`, `prop_top_10_percent`, `prop_top_25_percent`, `prop_top_50_percent`, `coeff_of_var_ratio` |
| `eval_expr` | `med_prop_na`, `med_prop_above_med` |
| `eval_compactness` | `autocor_median` |
| `eval_stan` | `standardization_comp` |
| `compare_metrics` | `rho_mean_med`, `rho_pca1_med`, `rho_mean_pca1`, `prop_pca1_var` |

### Step 2: Isolate the Module

Run just the failing module on the same input:

```python
# Python
from pysigqc.eval_var import compute_var
result = compute_var(signatures, names_sigs, datasets, names_datasets)
print(result["radar_values"]["sig_name"]["dataset_name"])
```

```r
# R (source R_refactored/)
for (f in list.files("R_refactored", pattern="[.]R$", full.names=TRUE)) source(f)
result <- compute_var(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)
result$radar_values[["sig_name"]][["dataset_name"]]
```

### Step 3: Check Common Causes

1. **NA handling:** R's `na.omit()` on matrices drops entire rows; Python's `np.nanmean()` skips individual NAs. The Python port uses `dropna(axis=0, how="any")` to match R.

2. **Sort behavior:** R's `sort()` removes NAs by default (`na.last=NA`). Python's `np.sort()` puts NaN at end. The Python port uses `np.nanmedian()` where R uses `median(sort(...))`.

3. **Correlation with constant columns:** R's `cor()` returns NaN per-pair for constant columns. `scipy.stats.spearmanr` on a full matrix may return scalar NaN. The Python port uses pairwise loops.

4. **PCA sign:** R's `prcomp` and sklearn's `PCA` may flip eigenvector signs. Compare `abs()` values for PCA-related metrics.

5. **Gaussian mixture:** R's `mclust::Mclust()` and sklearn's `GaussianMixture` use different EM initializations. Allow wider tolerance (~5%).

### Step 4: Trace the Intermediate Values

Both R and Python compute functions return intermediate data (not just radar values). Compare these:

- `eval_var`: `mean_sd_tables`, `cv_values`, `skewness_values`
- `eval_expr`: `na_proportions`, `prop_expressed`
- `eval_compactness`: `autocor_matrices`, `rank_products`
- `compare_metrics`: `scores`, `correlations`, `pca_results`
- `eval_stan`: `raw_vs_z_cors`

---

## Python Package Usage

### Basic Usage

```python
import pandas as pd
from pysigqc.eval_var import compute_var
from pysigqc.eval_expr import compute_expr
from pysigqc.eval_compactness import compute_compactness
from pysigqc.eval_stan import compute_stan
from pysigqc.compare_metrics import compute_metrics
from pysigqc.radar_chart import compute_radar

# Load data
dataset_A = pd.read_csv("expression_A.csv", index_col=0)  # genes x samples
datasets = {"dataset_A": dataset_A}
names_datasets = ["dataset_A"]

# Define signatures
signatures = {
    "my_sig": ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"]
}
names_sigs = ["my_sig"]

# Run individual modules
var_result = compute_var(signatures, names_sigs, datasets, names_datasets)
expr_result = compute_expr(signatures, names_sigs, datasets, names_datasets)
compact_result = compute_compactness(signatures, names_sigs, datasets, names_datasets)
stan_result = compute_stan(signatures, names_sigs, datasets, names_datasets)
metrics_result = compute_metrics(signatures, names_sigs, datasets, names_datasets)

# Assemble radar chart
radar_values = {}
for sig in names_sigs:
    radar_values[sig] = {}
    for ds in names_datasets:
        vals = {}
        vals.update(var_result["radar_values"][sig][ds])
        vals.update(expr_result["radar_values"][sig][ds])
        vals.update(compact_result["radar_values"][sig][ds])
        vals.update(metrics_result["radar_values"][sig][ds])
        vals.update(stan_result["radar_values"][sig][ds])
        radar_values[sig][ds] = vals

radar = compute_radar(radar_values, names_sigs, names_datasets)
print(radar["output_table"])  # DataFrame with 14 metrics per sig/dataset combo
```

### Negative Control Analysis

```python
from pysigqc.negative_control import run_negative_control

result = run_negative_control(
    signatures={"my_sig": ["GENE1", "GENE2", "GENE3"]},
    datasets={"dataset_A": dataset_A},
    num_resampling=100,
    seed=42,
)

# result["negative_controls"]["dataset_A"]["my_sig"]["summary"]
#   → DataFrame with mean + quantile rows for 14 metrics
# result["permutation_controls"]["dataset_A"]["my_sig"]["summary"]
#   → Same for permutation-based controls
```

### Input Requirements

- **Expression matrices:** pandas DataFrames with genes as rows (index) and samples as columns. Should be normalized, batch-corrected, and log-transformed.
- **Signatures:** dict mapping signature names to lists of gene names (must match expression matrix row index).
- **Minimum:** 2 genes per signature, 2 samples per dataset.

---

## Key R/Python Behavioral Differences

These are the specific R behaviors that required careful matching in the Python port. Relevant if debugging discrepancies or extending the code.

| R Behavior | Python Equivalent | Where It Matters |
|-----------|-------------------|-----------------|
| `na.omit(matrix)` drops entire rows | `df.dropna(axis=0, how="any")` | `eval_expr` threshold computation |
| `sort(x)` removes NAs (`na.last=NA`) | `np.sort(x[~np.isnan(x)])` | `eval_expr` proportion calculations |
| `cor(matrix)` returns NaN per-pair for constant cols | Pairwise `spearmanr` loop | `eval_compactness` autocorrelation |
| `prcomp()` eigenvector sign | sklearn `PCA` may flip sign | `compare_metrics` PCA correlations |
| `moments::skewness()` uses biased formula | `scipy.stats.skew(bias=True)` | `eval_var` skewness ratio |
| `Mclust()` BIC-based model selection | `GaussianMixture` with manual BIC | `compare_metrics` mixture models |
| `sd()` uses `n-1` denominator | `np.std(ddof=1)` | All z-transforms |

---

## Test Coverage Summary

### R Tests (testthat)

9 test files covering:
- `test-compute_var.R` — structure, 6 radar metrics, bounds
- `test-compute_expr.R` — structure, NA detection, expression bounds
- `test-compute_compactness.R` — autocorrelation matrix properties, diagonal=1
- `test-compute_stan.R` — correlation bounds, zero-variance handling
- `test-compute_metrics.R` — PCA variance, score lengths, mixture models
- `test-compute_struct.R` — NA padding (BUG-8), biclust, z-scores
- `test-compute_radar.R` — output table shape, metric assembly
- `test-negative_control.R` — QC metrics structure, file I/O
- `test-integration.R` — full pipeline, 14 metrics, determinism

### Python Tests (pytest)

9 test files:
- `test_eval_var.py` — 7 tests + cross-validation against R references
- `test_eval_expr.py` — 7 tests + cross-validation
- `test_eval_compactness.py` — 6 tests + cross-validation
- `test_eval_stan.py` — 5 tests + cross-validation
- `test_compare_metrics.py` — 7 tests + cross-validation
- `test_eval_struct.py` — 6 tests
- `test_radar_chart.py` — 6 tests + cross-validation
- `test_negative_control.py` — 4 tests (structure, determinism)
- `test_medium_cross_validation.py` — 6 tests on medium fixture (100 genes, 50 samples)
