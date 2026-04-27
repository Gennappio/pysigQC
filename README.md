# pysigQC

Python port of sigQC, an R package for quality control of gene signatures across genomic datasets. Based on the methodology published in Nature Protocols (Dhawan et al., 2019).

## Quick Start

```python
from pysigqc_joblib.pipeline import run_pipeline
from pysigqc.plots import plot_all

# Prepare your data
signatures = {
    "my_signature": ["GENE1", "GENE2", "GENE3", ...],
}
datasets = {
    "dataset_A": expression_matrix_A,  # pandas DataFrame, genes x samples
    "dataset_B": expression_matrix_B,
}

# Run the QC pipeline
result = run_pipeline(
    signatures, 
    list(signatures.keys()), 
    datasets, 
    list(datasets.keys())
)

# Generate all plots
plot_paths = plot_all(result, out_dir="qc_output/")
```

## Documentation Files

This repository contains several documentation files. Here's what each one covers:

| File | Description |
|------|-------------|
| **README.md** | This file. Overview and quick start guide |
| **[PLOTS_GUIDE.md](PLOTS_GUIDE.md)** | Detailed explanation of each output plot and how to interpret it |
| **[GUIDE.md](GUIDE.md)** | Comprehensive technical guide covering architecture, API reference, and implementation details |
| **[sigQC-master/AUDIT.md](sigQC-master/AUDIT.md)** | Code audit documenting bugs, issues, and corrections in the original R code |
| **[sigQC-master/README.md](sigQC-master/README.md)** | Original R package documentation |

## Project Structure

```
pysigQC/
├── pysigqc/                 # Core Python implementation
│   ├── __init__.py
│   ├── variability.py       # Mean/SD, CV, skewness analysis
│   ├── expression.py        # NA proportions, expression levels
│   ├── compactness.py       # Autocorrelation analysis
│   ├── standardization.py   # Raw vs z-transformed comparison
│   ├── structure.py         # Clustering, biclustering
│   ├── metrics.py           # Scoring methods (mean, median, PCA1)
│   ├── radar_chart.py       # Summary radar chart computation
│   ├── negative_control.py  # Permutation controls
│   └── plots/               # Plotting layer (decoupled)
│       ├── plot_var.py
│       ├── plot_expr.py
│       ├── plot_compactness.py
│       ├── plot_struct.py
│       ├── plot_metrics.py
│       ├── plot_stan.py
│       ├── plot_radar.py
│       └── plot_negative_control.py
│
├── pysigqc_joblib/          # Parallelized version using joblib
│   ├── pipeline.py          # Main entry point
│   └── ...
│
├── sigQC-master/            # Original R package (reference)
│   ├── R/                   # Original R source
│   ├── R_corrected/         # Bug-fixed R source
│   └── AUDIT.md             # Detailed bug documentation
│
├── PLOTS_GUIDE.md           # Plot interpretation guide
├── GUIDE.md                 # Technical documentation
└── README.md                # This file
```

## Output Files

Running the pipeline generates these PDF files in your output directory:

| Category | File | Description |
|----------|------|-------------|
| **Summary** | `sig_radarplot.pdf` | 14-metric radar chart overview |
| **Variability** | `sig_mean_vs_sd.pdf` | Mean vs SD scatter plots |
| **Expression** | `sig_expr_barcharts.pdf` | Expression proportion bars |
| | `sig_expr_barcharts_NA_values.pdf` | Missing value bars |
| | `sig_expr_density_plots.pdf` | Expression density curves |
| **Compactness** | `sig_autocor_hmps.pdf` | Gene-gene correlation heatmaps |
| | `sig_autocor_dens.pdf` | Autocorrelation density overlay |
| **Structure** | `sig_heatmaps_hmps.pdf` | Clustered expression heatmaps |
| | `sig_biclust_*.pdf` | Biclustering results |
| **Metrics** | `sig_compare_metrics_<sig>.pdf` | Scoring method comparison |
| | `sig_qq_plots_<sig>.pdf` | QQ normality plots |
| **Standardization** | `sig_standardisation_comp.pdf` | Raw vs z-score comparison |

See **[PLOTS_GUIDE.md](PLOTS_GUIDE.md)** for detailed interpretation of each plot.

## Key Improvements Over R Version

1. **Decoupled computation and plotting** — Run computation once, generate plots separately
2. **Parallel processing** — Uses joblib for multi-core computation
3. **Bug fixes** — Several critical bugs in R code have been fixed (see AUDIT.md)
4. **Better error handling** — Explicit error handling instead of silent failures

## Requirements

- Python 3.10+
- numpy, pandas, scipy, scikit-learn
- matplotlib, seaborn
- joblib (for parallelization)

## Installing the R reference env (optional)

The Python port runs standalone. The R reference stack is only needed for the cross-validation smoke test (`TEST_stud/test3.py`) and for running the original sigQC against your own data. It is fully described by [`environment-r.yml`](environment-r.yml) and installed via conda:

```bash
./install_r.sh                 # creates conda env "sigqc-r" from environment-r.yml
./install_r.sh my-env-name     # custom env name
./install_r.sh --update        # refresh an existing env
```

R itself, every CRAN dependency, and every Bioconductor dependency (GSVA, ComplexHeatmap, biomaRt, RankProd, ...) come from `conda-forge` + `bioconda` — no system libraries or `BiocManager` calls required. See [TESTS_PYTHON.md](TESTS_PYTHON.md) for how to point `test3.py` at the resulting env without activating it.

## Citation

If you use sigQC in your research, please cite:

> Dhawan A, et al. "Guidelines for using sigQC for systematic evaluation of gene signatures." Nature Protocols (2019).

## License

See the original sigQC package for license terms.
