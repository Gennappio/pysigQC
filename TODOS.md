# TODOS — Testing, Discrepancy Analysis, and Performance Comparison

Actionable checklist for validating all code versions, investigating R/Python discrepancies, and benchmarking `pysigqc` vs `pysigqc_joblib` at increasing scale.

---

## 1. Generate Test Fixtures (R side)

All subsequent steps depend on these fixtures existing.

- [ ] **1.1** Run small fixture generator on the cluster:
  ```bash
  cd sigQC-master
  Rscript tests/fixtures/fixture_generator.R
  ```
  Produces: `fixture_small.rds`, `fixture_dataset_A.csv`, `fixture_dataset_B.csv`, `fixture_signatures.csv`

- [ ] **1.2** Run medium fixture generator:
  ```bash
  Rscript tests/fixtures/fixture_medium_generator.R
  ```
  Produces: `fixture_medium.rds`, `fixture_medium_dataset_{A,B,C}.csv`, `fixture_medium_signatures.csv`

- [ ] **1.3** Generate small reference outputs:
  ```bash
  Rscript tests/fixtures/generate_reference_outputs.R
  ```
  Produces: 43 CSV files in `tests/fixtures/reference_outputs/`

- [ ] **1.4** Generate medium reference outputs:
  ```bash
  Rscript tests/fixtures/generate_medium_reference.R
  ```
  Produces: 76 CSV files in `tests/fixtures/reference_outputs_medium/`

- [ ] **1.5** Verify file counts:
  ```bash
  ls tests/fixtures/reference_outputs/ | wc -l    # expect 43
  ls tests/fixtures/reference_outputs_medium/ | wc -l  # expect 76
  ```

**SLURM shortcut:** `sbatch slurm/run_r_tests.slurm` runs steps 1.1–1.4 plus R unit tests in one job.

---

## 2. Run R Unit Tests (testthat)

- [ ] **2.1** Run the full R test suite:
  ```bash
  cd sigQC-master
  Rscript tests/testthat.R
  ```

- [ ] **2.2** Verify all 9 test files pass:
  - `test-compute_var.R`
  - `test-compute_expr.R`
  - `test-compute_compactness.R`
  - `test-compute_stan.R`
  - `test-compute_metrics.R`
  - `test-compute_struct.R`
  - `test-compute_radar.R`
  - `test-negative_control.R`
  - `test-integration.R`

- [ ] **2.3** If any test fails, record:
  - Which test file and test name
  - The assertion that failed
  - The actual vs expected values
  - Whether the failure is in the `R_refactored/` code or the test expectations

---

## 3. Run Python Unit Tests (pysigqc)

- [ ] **3.1** Install the package:
  ```bash
  pip install -e ".[dev]"
  ```

- [ ] **3.2** Run all Python tests:
  ```bash
  pytest tests/ -v --tb=short
  ```

- [ ] **3.3** Verify all 65 tests pass (9 test files + 1 equivalence file):
  - `test_eval_var.py` — 7 tests
  - `test_eval_expr.py` — 7 tests
  - `test_eval_compactness.py` — 6 tests
  - `test_eval_stan.py` — 5 tests
  - `test_compare_metrics.py` — 7 tests
  - `test_eval_struct.py` — 6 tests
  - `test_radar_chart.py` — 6 tests
  - `test_negative_control.py` — 4 tests
  - `test_medium_cross_validation.py` — 6 tests
  - `test_joblib_equivalence.py` — 11 tests

- [ ] **3.4** Run only cross-validation tests (R reference required):
  ```bash
  pytest tests/ -v -k "cross_validation or medium"
  ```

- [ ] **3.5** Run only unit tests (no R reference needed):
  ```bash
  pytest tests/ -v -k "not cross_validation and not medium"
  ```

**SLURM shortcut:** `sbatch slurm/run_py_tests.slurm`

---

## 4. Investigate R vs Python Discrepancies

For each metric that fails cross-validation or shows unexpected values.

### 4.1 Per-Module Discrepancy Checklist

For each module, compare R and Python outputs on the **same** fixture data:

- [ ] **4.1.1 eval_var** — 6 metrics:
  ```python
  from pysigqc.eval_var import compute_var
  result = compute_var(signatures, names_sigs, datasets, names_datasets)
  for sig in names_sigs:
      for ds in names_datasets:
          print(sig, ds, result["radar_values"][sig][ds])
  ```
  Compare against `reference_outputs/var_radar_{sig}_{ds}.csv`.
  Known pitfalls: `scipy.stats.skew(bias=True)` must match R's `moments::skewness`.

- [ ] **4.1.2 eval_expr** — 2 metrics:
  Check threshold computation: Python must use `dropna(axis=0, how="any")` to match R's `na.omit()` on matrices (drops entire rows, not individual NAs).
  Check `med_prop_above_med`: Python must set `gene_expr_props[has_na] = np.nan` to match R's NA propagation from `rowSums`.

- [ ] **4.1.3 eval_compactness** — 1 metric:
  Check that pairwise `spearmanr` is used (not full-matrix `spearmanr`), because R's `cor()` handles constant columns per-pair. Verify constant-gene pairs produce NaN in both.

- [ ] **4.1.4 eval_stan** — 1 metric:
  Check zero-variance genes: Python's `z_transform` should return zeros (not NaN/Inf). Verify `np.nanstd(ddof=1)` matches R's `sd()`.

- [ ] **4.1.5 compare_metrics** — 4 metrics:
  PCA sign ambiguity: compare `abs()` values for `rho_pca1_med` and `rho_mean_pca1`.
  GaussianMixture vs Mclust: allow 5% tolerance. Different EM initializations produce different cluster assignments.

- [ ] **4.1.6 Full radar** — 14 metrics:
  Run `test_medium_cross_validation.py::test_full_radar_medium` and inspect any failures.

### 4.2 Systematic Discrepancy Trace

When a specific metric disagrees beyond tolerance:

- [ ] **4.2.1** Identify the exact (signature, dataset, metric) triple.
- [ ] **4.2.2** Print the R reference value from the CSV file.
- [ ] **4.2.3** Print the Python computed value.
- [ ] **4.2.4** Check the intermediate values (not just the final radar metric):
  ```python
  # Example for eval_var
  result = compute_var(sigs, names_sigs, datasets, names_datasets)
  print("SD table:", result["mean_sd_tables"][sig][ds])
  print("All SD:", result["all_sd"][sig][ds].head(20))
  ```
- [ ] **4.2.5** Run the same computation in R interactively:
  ```r
  source("R_refactored/eval_var_loc.R")
  result <- compute_var(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)
  result$radar_values[["sig_name"]][["dataset_name"]]
  ```
- [ ] **4.2.6** Compare intermediate values step by step until the divergence point is found.
- [ ] **4.2.7** Document the root cause (see known causes in GUIDE.md section "Key R/Python Behavioral Differences").
- [ ] **4.2.8** If a genuine Python bug: fix in `pysigqc/`, re-run cross-validation, then propagate the fix to `pysigqc_joblib/`.

---

## 5. Verify pysigqc vs pysigqc_joblib Numerical Equivalence

- [ ] **5.1** Run the dedicated equivalence tests:
  ```bash
  pytest tests/test_joblib_equivalence.py -v
  ```

- [ ] **5.2** Verify all 11 tests pass:
  - Small fixture: var, expr, compactness, stan, metrics, full_radar (6 tests)
  - Medium fixture: var, expr, compactness, stan, metrics (5 tests)

- [ ] **5.3** If a test fails, the discrepancy is a **bug in pysigqc_joblib** (the reference pysigqc is ground truth). Investigate:
  - Which module?
  - Which metric?
  - Is it a constant-gene edge case (the `fill_diagonal` fix)?
  - Is it a NaN handling difference?
  - Is it a floating-point order-of-operations difference (acceptable if < 1e-10)?

- [ ] **5.4** Test with parallel execution:
  ```python
  from pysigqc_joblib.pipeline import run_pipeline
  r1 = run_pipeline(sigs, names_sigs, datasets, names_datasets, n_jobs=1)
  r2 = run_pipeline(sigs, names_sigs, datasets, names_datasets, n_jobs=-1)
  # Compare r1 and r2 radar tables — should be identical
  ```

---

## 6. Performance Comparison: pysigqc vs pysigqc_joblib

### 6.1 Medium Fixture Benchmark (100 genes, 50 samples)

- [ ] **6.1.1** Run the per-module timing comparison:
  ```python
  import pandas as pd
  from pathlib import Path

  FIXTURES_DIR = Path("sigQC-master/tests/fixtures")
  datasets = {name: pd.read_csv(FIXTURES_DIR / f"fixture_medium_{name}.csv", index_col=0)
               for name in ["dataset_A", "dataset_B", "dataset_C"]}
  df = pd.read_csv(FIXTURES_DIR / "fixture_medium_signatures.csv")
  sigs = {sig: group["gene"].tolist() for sig, group in df.groupby("signature")}
  names_sigs = list(sigs.keys())
  names_datasets = list(datasets.keys())

  import pysigqc.eval_var as ref_var
  import pysigqc_joblib.eval_var as opt_var
  # ... repeat for each module

  ref = ref_var.compute_var(sigs, names_sigs, datasets, names_datasets)
  opt = opt_var.compute_var(sigs, names_sigs, datasets, names_datasets)
  print(f"pysigqc: {ref['elapsed_seconds']:.4f}s")
  print(f"joblib:  {opt['elapsed_seconds']:.4f}s")
  print(f"speedup: {ref['elapsed_seconds'] / opt['elapsed_seconds']:.1f}x")
  ```

- [ ] **6.1.2** Record results for all modules:

  | Module | pysigqc (s) | pysigqc_joblib (s) | Speedup |
  |--------|-------------|-------------------|---------|
  | eval_var | | | |
  | eval_expr | | | |
  | eval_compactness | | | |
  | eval_stan | | | |
  | compare_metrics | | | |
  | radar_chart | | | |

### 6.2 Synthetic Large-Scale Benchmark

Create progressively larger synthetic datasets to measure scaling behavior.

- [ ] **6.2.1** Generate synthetic datasets at increasing scale:
  ```python
  import numpy as np
  import pandas as pd

  def make_synthetic(n_genes, n_samples, seed=42):
      rng = np.random.default_rng(seed)
      arr = rng.normal(5, 2, size=(n_genes, n_samples))
      genes = [f"GENE{i:05d}" for i in range(n_genes)]
      samples = [f"S{i:06d}" for i in range(n_samples)]
      return pd.DataFrame(arr, index=genes, columns=samples)

  # Signature: 20 genes from the first 20 rows
  sigs = {"test_sig": [f"GENE{i:05d}" for i in range(20)]}
  names_sigs = ["test_sig"]
  ```

- [ ] **6.2.2** Benchmark at each scale (single dataset):

  | Scale | n_genes | n_samples | pysigqc (s) | joblib (s) | Speedup | Peak RAM |
  |-------|---------|-----------|-------------|-----------|---------|----------|
  | Small | 100 | 50 | | | | |
  | Medium | 1,000 | 1,000 | | | | |
  | Large | 10,000 | 10,000 | | | | |
  | XL | 20,000 | 100,000 | | | | |
  | SC-scale | 20,000 | 1,000,000 | | | | |

- [ ] **6.2.3** For each scale, measure per-module and total pipeline time:
  ```python
  from pysigqc_joblib.pipeline import run_pipeline
  result = run_pipeline(sigs, names_sigs, datasets, names_datasets, n_jobs=-1)
  print(f"Total: {result['elapsed_seconds']:.2f}s")
  print(f"  var:        {result['var_result']['elapsed_seconds']:.2f}s")
  print(f"  expr:       {result['expr_result']['elapsed_seconds']:.2f}s")
  print(f"  compactness:{result['compact_result']['elapsed_seconds']:.2f}s")
  print(f"  stan:       {result['stan_result']['elapsed_seconds']:.2f}s")
  print(f"  metrics:    {result['metrics_result']['elapsed_seconds']:.2f}s")
  ```

- [ ] **6.2.4** Measure peak memory at XL and SC-scale:
  ```python
  import tracemalloc
  tracemalloc.start()
  result = run_pipeline(sigs, names_sigs, datasets, names_datasets, n_jobs=1)
  current, peak = tracemalloc.get_traced_memory()
  tracemalloc.stop()
  print(f"Peak memory: {peak / 1024**2:.1f} MB")
  ```

### 6.3 Negative Control Benchmark

This is the most important benchmark because negative controls dominate total runtime.

- [ ] **6.3.1** Compare sequential negative control timing:
  ```python
  from pysigqc.negative_control import run_negative_control as nc_ref
  from pysigqc_joblib.negative_control import run_negative_control as nc_opt

  # num_resampling=50, n_jobs=1 (sequential for fair comparison)
  ref = nc_ref(sigs_1sig, datasets_1ds, num_resampling=50, seed=42)
  opt = nc_opt(sigs_1sig, datasets_1ds, num_resampling=50, n_jobs=1, seed=42)
  print(f"Ref: {ref['elapsed_seconds']:.2f}s, Opt: {opt['elapsed_seconds']:.2f}s")
  ```

- [ ] **6.3.2** Measure parallel speedup vs core count:

  | n_jobs | Time (s) | Speedup vs n_jobs=1 |
  |--------|----------|---------------------|
  | 1 | | 1.0x |
  | 2 | | |
  | 4 | | |
  | 8 | | |
  | -1 (all) | | |

- [ ] **6.3.3** Test at scale: num_resampling=1000 on 10k-sample dataset:
  ```python
  ds_10k = {"ds": make_synthetic(5000, 10000)}
  result = nc_opt(sigs, ds_10k, num_resampling=1000, n_jobs=-1, seed=42)
  print(f"1000 resamplings on 10k samples: {result['elapsed_seconds']:.1f}s")
  ```

---

## 7. BRCA Paper Data Cross-Validation

Validate on the actual data from the Nature Protocols paper.

- [ ] **7.1** Ensure BRCA data is placed in `sigQC-master/input_files/`:
  - Expression matrices as `exprs*.csv`
  - Signatures as `sig*.csv` (long format: signature, gene columns)

- [ ] **7.2** Run the BRCA cross-validation SLURM job:
  ```bash
  cd sigQC-master/slurm
  sbatch run_brca_cross_validation.slurm
  ```

- [ ] **7.3** Check the output:
  ```bash
  cat brca_cross_validation/py_output/brca_mismatches.csv
  ```

- [ ] **7.4** For each mismatch reported:
  - [ ] Is it a PCA sign issue? (check if `pca1` is in the column name)
  - [ ] Is it a GaussianMixture tolerance issue? (check if diff < 5%)
  - [ ] Is it a genuine discrepancy? (follow procedure in section 4.2)

- [ ] **7.5** Compare R and Python radar output tables visually:
  ```python
  import pandas as pd
  r_table = pd.read_csv("brca_cross_validation/r_output/radar_output_table.csv", index_col=0)
  py_table = pd.read_csv("brca_cross_validation/py_output/radar_output_table.csv", index_col=0)
  diff = (py_table - r_table).abs()
  print("Max absolute difference per metric:")
  print(diff.max())
  ```

- [ ] **7.6** Run BRCA data through `pysigqc_joblib` and compare with `pysigqc`:
  ```python
  from pysigqc_joblib.pipeline import run_pipeline
  result = run_pipeline(sigs, names_sigs, datasets, names_datasets, n_jobs=-1)
  # Compare result["radar_result"]["output_table"] with pysigqc output
  ```

---

## 8. Edge Case Stress Tests

Verify robustness on pathological inputs.

- [ ] **8.1** All-NA signature gene (every sample is NaN for one gene):
  ```python
  ds = datasets["dataset_A"].copy()
  ds.loc["GENE001"] = np.nan  # make one sig gene all-NaN
  # Run each module, check no crash, NaN propagates correctly
  ```

- [ ] **8.2** Single-gene signature (k=1, minimum viable):
  ```python
  sigs_1gene = {"tiny_sig": ["GENE001"]}
  # Compactness should return autocor_median=0.0 (no pairs)
  # PCA should be skipped (need >=2 genes)
  ```

- [ ] **8.3** All-constant expression matrix (zero variance everywhere):
  ```python
  ds_const = pd.DataFrame(5.0, index=genes, columns=samples)
  # z_transform should return zeros, no Inf/NaN
  # Spearman correlations should be NaN
  ```

- [ ] **8.4** Signature with no genes present in dataset:
  ```python
  sigs_missing = {"ghost_sig": ["NONEXISTENT_1", "NONEXISTENT_2", "NONEXISTENT_3"]}
  # gene_intersection returns [], modules should handle gracefully
  ```

- [ ] **8.5** Very large signature (k=500 genes):
  ```python
  sigs_large = {"big_sig": [f"GENE{i:05d}" for i in range(500)]}
  # Compactness: 500x500 correlation matrix = 250k pairs
  # Should still complete in reasonable time with corrcoef optimization
  ```

- [ ] **8.6** Dataset with duplicate gene names:
  ```python
  # pandas allows duplicate index — check modules don't produce wrong results
  ds_dup = pd.concat([datasets["dataset_A"], datasets["dataset_A"].iloc[:5]])
  ```

---

## 9. Reproducibility Checks

- [ ] **9.1** Determinism of pysigqc (same inputs, same outputs, multiple runs):
  ```python
  r1 = compute_var(sigs, names_sigs, datasets, names_datasets)
  r2 = compute_var(sigs, names_sigs, datasets, names_datasets)
  assert r1["radar_values"] == r2["radar_values"]
  ```

- [ ] **9.2** Determinism of pysigqc_joblib with n_jobs=1:
  ```python
  r1 = nc_opt(sigs, datasets, num_resampling=10, n_jobs=1, seed=42)
  r2 = nc_opt(sigs, datasets, num_resampling=10, n_jobs=1, seed=42)
  pd.testing.assert_frame_equal(
      r1["negative_controls"]["dataset_A"]["sig"]["metrics_table"],
      r2["negative_controls"]["dataset_A"]["sig"]["metrics_table"],
  )
  ```

- [ ] **9.3** Determinism of pysigqc_joblib with n_jobs=-1 (parallel):
  ```python
  # Same seed should give same results even with parallel execution
  r1 = nc_opt(sigs, datasets, num_resampling=10, n_jobs=-1, seed=42)
  r2 = nc_opt(sigs, datasets, num_resampling=10, n_jobs=-1, seed=42)
  # Note: parallel RNG uses SeedSequence.spawn(), so results will differ
  # from n_jobs=1, but should be self-consistent across multiple runs
  ```

- [ ] **9.4** Cross-platform reproducibility:
  Run the same fixture through both Python packages on:
  - [ ] Local machine (macOS)
  - [ ] SLURM cluster (Linux)
  Verify radar values match within floating-point tolerance (rtol=1e-12).

---

## 10. Documentation and Cleanup

- [ ] **10.1** After all tests pass, update GUIDE.md with:
  - Final test count
  - Any new edge cases discovered
  - Benchmark numbers from the cluster

- [ ] **10.2** Update OPTIMIZATIONS.md with:
  - Actual benchmark numbers at each scale
  - Memory measurements
  - Any optimizations added during testing

- [ ] **10.3** Verify all files are tracked in git:
  ```bash
  git status
  ```

- [ ] **10.4** Review for any hardcoded paths or credentials before committing.

- [ ] **10.5** Final full test run:
  ```bash
  pytest tests/ -v --tb=short
  ```
  All 65 tests must pass.
