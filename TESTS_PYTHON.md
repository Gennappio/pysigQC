# Python Test Suite — Explained

This document explains every Python test file, what question each test is asking, and what a failure means. Tests live in `tests/`. Run with:

```bash
pytest tests/ -v --tb=short
```

For fixture and output file definitions, see `FIXTURES.md`.
For the metrics being tested and their biological meaning, see `FIXTURES.md § Module prefixes`.
For the bugs these tests guard against, see `sigQC-master/AUDIT.md`.

---

## conftest.py — Shared fixtures

Not a test file, but defines the pytest fixtures (small dataset) used by all tests:

- `signatures` — dict mapping signature name → list of gene IDs
- `datasets` — dict mapping dataset name → pandas DataFrame (genes × samples)
- `names_sigs`, `names_datasets` — ordered name lists
- `ref_dir` — path to `sigQC-master/tests/fixtures/reference_outputs/`

These fixtures are loaded from `fixture_dataset_A.csv`, `fixture_dataset_B.csv`, and `fixture_signatures.csv`. If those files don't exist, all cross-validation tests are skipped (not failed).

---

## test_eval_var.py — Variability module

Tests `pysigqc.eval_var.compute_var()`, which measures how variable signature genes are relative to the full transcriptome.

| Test | What it checks | Why it matters |
|------|---------------|----------------|
| `test_compute_var_structure` | Return dict has exactly the expected keys | Guards against accidentally renaming or removing outputs |
| `test_compute_var_6_metrics` | The radar_values sub-dict has exactly 6 metric keys | Catches metric additions/removals that would silently change the radar chart |
| `test_compute_var_values_in_0_1` | All 6 metric values are in [0, 1] | The radar chart normalises everything to [0,1]; values outside this range corrupt the visualisation |
| `test_compute_var_mean_sd_table_dims` | The mean/SD table has one row per gene in the intersection, two columns (Mean, SD) | The intersection is genes present in both the signature and the dataset; `gene_missing` must not appear |
| `test_compute_var_missing_gene_excluded` | `gene_missing` is absent from `inter["compact_sig"]["dataset_A"]`, leaving exactly 4 genes | Directly tests the missing-gene handling: a gene in the signature but not the dataset must be silently dropped, not raise an error |
| `test_compute_var_deterministic` | Running twice with the same inputs gives identical outputs | `eval_var` has no randomness; any non-determinism would indicate a hidden global state bug |
| `test_compute_var_cross_validation` | Every metric value matches the R reference within `rtol=1e-4` | The primary correctness test; verifies the Python port computes the same numbers as the R original |

---

## test_eval_expr.py — Expression module

Tests `pysigqc.eval_expr.compute_expr()`, which measures how well expressed signature genes are in each dataset.

| Test | What it checks | Why it matters |
|------|---------------|----------------|
| `test_structure` | Return keys | Structure guard |
| `test_2_metrics` | Exactly `med_prop_na` and `med_prop_above_med` present | These are the two expression radar metrics |
| `test_values_in_0_1` | Both proportions in [0, 1] | They are proportions; values outside [0,1] indicate a division error |
| `test_detects_na` | `gene_missing` has `na_proportion == 1.0` in `dataset_A` | A gene absent from the dataset should have 100% NA, not 0% or a crash |
| `test_thresholds_positive` | The auto-computed expression threshold is > 0 | The threshold is the median expression across all genes; it should always be positive for log-expression data |
| `test_custom_thresholds` | If thresholds are supplied manually, they override auto-computation | This is the user-facing API for per-dataset thresholds |
| `test_cross_validation` | Matches R reference within `rtol=1e-4` | Correctness test; key pitfall: R's `na.omit()` on matrices drops entire rows, so Python must use `dropna(axis=0, how="any")` not per-cell NA removal |

---

## test_eval_compactness.py — Compactness module

Tests `pysigqc.eval_compactness.compute_compactness()`, which measures how co-expressed signature genes are with each other (pairwise Spearman correlation).

| Test | What it checks | Why it matters |
|------|---------------|----------------|
| `test_structure` | Return keys | Structure guard |
| `test_autocor_median_present` | The radar metric key `autocor_median` exists | This is the single most important compactness metric |
| `test_compact_higher_than_diffuse` | `compact_sig.autocor_median > diffuse_sig.autocor_median` for every dataset | The small fixture was designed so this must hold: `compact_sig` genes are generated from a common signal, `diffuse_sig` genes are independent. If this fails, the correlation computation is wrong. |
| `test_autocor_matrix_symmetric` | The pairwise correlation matrix is symmetric and has 1.0 on the diagonal | Spearman correlation matrices are always symmetric; diagonal must be 1 (a gene is perfectly correlated with itself). This guards the `fill_diagonal` fix — `np.corrcoef` returns NaN on the diagonal for constant-variance rows, which would corrupt `autocor_median`. |
| `test_autocor_values_bounded` | All finite matrix values are in [-1, 1] | Correlation is bounded; values outside this range indicate numerical instability |
| `test_cross_validation` | Matches R reference | Correctness test; the key technical detail is that R's `cor()` handles constant columns per-pair whereas numpy's `corrcoef` does not — see `OPTIMIZATIONS.md § spearman_matrix` |

---

## test_eval_stan.py — Standardisation module

Tests `pysigqc.eval_stan.compute_stan()`, which asks: does z-transforming each gene's expression change the sample ranking produced by the signature?

| Test | What it checks | Why it matters |
|------|---------------|----------------|
| `test_structure` | Return keys | Structure guard |
| `test_metric_range` | `standardization_comp` is in [-1, 1] | It is a Spearman correlation; must be bounded |
| `test_score_lengths` | Both `med_scores` and `z_transf_scores` have 10 entries | The fixture has 10 samples; one score per sample |
| `test_zero_variance_gene` | Running `diffuse_sig` on `dataset_B` (which has constant `gene_5`) does not produce NaN | Directly tests the BUG-2 fix (`sigQC-master/AUDIT.md § BUG-2`). If the z-transform divides by zero, the score becomes NaN/Inf, then the Spearman correlation is NaN, and the test fails. |
| `test_cross_validation` | Matches R reference | Correctness test |

---

## test_compare_metrics.py — PCA/scoring module

Tests `pysigqc.compare_metrics.compute_metrics()`, which compares three ways to score samples (mean, median, PCA1) and checks how consistent they are.

| Test | What it checks | Why it matters |
|------|---------------|----------------|
| `test_structure` | Return keys: `radar_values`, `scores`, `pca_results`, `score_cor_mats`, `mixture_models` | Structure guard |
| `test_4_radar_metrics` | Exactly the 4 PCA/correlation metrics present per (sig, ds) pair | These 4 metrics are separate from the 10 produced by var/expr/compact/stan |
| `test_correlation_values_bounded` | `rho_mean_med`, `rho_pca1_med`, `rho_mean_pca1` all in [-1, 1] | They are Spearman correlations |
| `test_prop_pca1_var_bounded` | `prop_pca1_var` in [0, 1] | It is a proportion of variance; must be in [0,1] |
| `test_score_lengths` | `med_scores` and `mean_scores` each have 10 entries (one per sample) | Sanity check on scoring output |
| `test_pca_variance_sums_to_1` | Sum of all principal component variances = 1.0 | Mathematical property of PCA; if this fails, the PCA computation has a normalisation bug |
| `test_cross_validation` | Matches R reference using `abs()` for PCA sign metrics | PCA eigenvectors have sign ambiguity: the first PC can point in either direction. `rho_pca1_med` may have opposite sign between R and Python but same absolute value. |

---

## test_eval_struct.py — Structure module

Tests `pysigqc.eval_struct.compute_struct()`, which performs biclustering and hierarchical clustering to look for block structure in gene expression.

| Test | What it checks | Why it matters |
|------|---------------|----------------|
| `test_structure` | Return keys | Structure guard |
| `test_all_row_names_union` | The union of gene names covers all genes present in any dataset | The struct module aligns matrices across datasets; every gene seen in at least one dataset should appear in the output |
| `test_padded_matrices_consistent_rows` | All dataset matrices have the same number of rows | Biclustering requires aligned matrices; row count must match across datasets |
| `test_padded_with_na_for_missing` | A gene present in dataset X but not dataset Y gets an all-NA row in Y | Tests the NA-padding fix (replaces the old BUG-8 min-value padding — see `sigQC-master/AUDIT.md § BUG-8`). An all-NA row is neutral; a min-value row introduces an artificial zero-variance gene. |
| `test_biclust_results_populated` | Each (sig, ds) has `z_scores`, `binarized`, `biclust_result`, `threshold` keys | The internal pipeline stages are all present |
| `test_any_biclusters_is_bool` | The `any_biclusters` flag is a Python bool | Used downstream to decide whether to draw bicluster plots |

---

## test_radar_chart.py — Radar chart module

Tests `pysigqc.radar_chart.compute_radar()`, which aggregates all 14 metrics into a summary table and computes radar chart areas.

| Test | What it checks | Why it matters |
|------|---------------|----------------|
| `test_structure` | Return keys: `radar_plot_mat`, `output_table`, `areas`, `legend_labels`, `radarplot_rownames` | Structure guard |
| `test_output_dimensions` | `output_table` has shape `(n_sigs × n_datasets, 14)` | 14 metrics, one row per (sig, ds) pair |
| `test_non_negative` | All values ≥ 0 | Radar chart values are normalised proportions; negatives would indicate a metric was not clipped before plotting |
| `test_areas_positive` | All radar chart areas ≥ 0 | Area is a sum of products; must be non-negative. Note: the area formula is a known approximation — see `sigQC-master/AUDIT.md § ISSUE-3`. |
| `test_fills_missing_metrics` | If only one metric is provided in `radar_values`, the output table still has 14 columns (zeros for the rest) | The radar chart always produces a fixed-width 14-column table regardless of which modules were run |
| `test_cross_validation` | Full radar table matches R reference, skipping the 4 `compare_metrics` columns if they are zero | Those 4 columns are zero when `compare_metrics` is excluded from the pipeline |

---

## test_negative_control.py — Negative control module

Tests `pysigqc.negative_control.run_negative_control()`, which estimates a null distribution for all 14 metrics by repeatedly drawing random gene sets.

| Test | What it checks | Why it matters |
|------|---------------|----------------|
| `test_compute_qc_metrics_structure` | `_compute_qc_metrics()` returns `radar_values` and `output_table` with correct shape | This is the internal helper that runs the full pipeline once; the NC module calls it many times |
| `test_compute_qc_metrics_non_negative` | Output table values ≥ 0 | Same as radar chart constraint |
| `test_run_negative_control_small` | With `num_resampling=3`: negative controls and permutation controls both have 3-row `metrics_table` and a 6-row summary (mean + 5 quantiles) | Verifies the resampling loop produces the right output shape. `num_resampling=3` is used to keep test time short. |
| `test_negative_control_deterministic` | Running with the same `seed=42` twice produces identical `metrics_table` | The BUG-1 fix (`sigQC-master/AUDIT.md § BUG-1`) replaced `runif()` with `sample()`; this test verifies the Python port uses proper reproducible random sampling. |

---

## test_medium_cross_validation.py — Medium fixture cross-validation

Runs each module on the 100 gene × 50 sample × 3 dataset × 5 signature medium fixture and compares to R reference outputs. This is the main integration test — if the small fixture tests pass but the medium tests fail, the bug is likely a numerical stability issue at realistic scale.

Each test (`test_var_medium`, `test_expr_medium`, etc.) runs one module on all 15 (sig, ds) pairs and checks every radar metric against the corresponding `reference_outputs_medium/` CSV.

`test_full_radar_medium` runs the entire pipeline end-to-end and compares the full 15-row × 14-column output table against `radar_output_table.csv`.

**These tests are skipped automatically if the medium fixture CSVs do not exist.** Generate them first — see `TODOS.md § 1`.

---

## test_joblib_equivalence.py — Optimised vs reference equivalence

Verifies that `pysigqc_joblib` (the parallel, vectorised version) produces **numerically identical results** to `pysigqc` (the reference version). See `OPTIMIZATIONS.md` for what was changed and why.

The file has two test classes:

**`TestSmallFixtureEquivalence`** — 6 tests on the small fixture:
- One test per module (`var`, `expr`, `compactness`, `stan`, `metrics`)
- One full pipeline test assembling all 14 metrics and comparing radar tables

**`TestMediumFixtureEquivalence`** — 5 tests on the medium fixture:
- One test per module (no full radar test at medium scale to keep runtime reasonable)

All tests use `np.testing.assert_allclose(rtol=1e-10)` — this is much tighter than the R cross-validation tolerance of `1e-4`. If this fails, `pysigqc_joblib` has a bug. The reference implementation `pysigqc` is always the ground truth.

The one exception is PCA metrics (`rho_pca1_*`, `rho_mean_pca1`), where `abs()` comparison is used because PCA sign ambiguity may differ between the two implementations even though the underlying computation is equivalent.
