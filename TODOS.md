# TODOS — Checklist

See supporting docs for context:
- `FIXTURES.md` — what the fixtures model and how to read output files
- `TESTS_PYTHON.md` — what each Python test checks and why
- `sigQC-master/AUDIT.md` — confirmed bugs and theoretical issues in the R code
- `OPTIMIZATIONS.md` — what was changed in `pysigqc_joblib` and why
- `GUIDE.md` — overall architecture, how to run things, R/Python differences

---

## 1. Generate fixtures (R side)

- [ ] **1.1** Small fixture:
  ```bash
  cd sigQC-master && Rscript tests/fixtures/fixture_generator.R
  ```
  Produces: `fixture_small.rds`, `fixture_dataset_A.csv`, `fixture_dataset_B.csv`, `fixture_signatures.csv`

- [ ] **1.2** Medium fixture:
  ```bash
  Rscript tests/fixtures/fixture_medium_generator.R
  ```
  Produces: `fixture_medium.rds`, `fixture_medium_dataset_{A,B,C}.csv`, `fixture_medium_signatures.csv`

- [ ] **1.3** Small reference outputs (R ground-truth values):
  ```bash
  Rscript tests/fixtures/generate_reference_outputs.R
  ```
  Produces: 43 CSVs in `tests/fixtures/reference_outputs/`

- [ ] **1.4** Medium reference outputs:
  ```bash
  Rscript tests/fixtures/generate_medium_reference.R
  ```
  Produces: 76 CSVs in `tests/fixtures/reference_outputs_medium/`

- [ ] **1.5** Verify counts:
  ```bash
  ls tests/fixtures/reference_outputs/ | wc -l       # expect 43
  ls tests/fixtures/reference_outputs_medium/ | wc -l # expect 76
  ```

**SLURM shortcut:** `sbatch slurm/run_r_tests.slurm` runs 1.1–1.4 plus the R unit tests.

---

## 2. Run R unit tests

- [ ] **2.1** Run full R test suite:
  ```bash
  cd sigQC-master && Rscript tests/testthat.R
  ```

- [ ] **2.2** All 9 test files must pass:
  `test-compute_var.R`, `test-compute_expr.R`, `test-compute_compactness.R`,
  `test-compute_stan.R`, `test-compute_metrics.R`, `test-compute_struct.R`,
  `test-compute_radar.R`, `test-negative_control.R`, `test-integration.R`

- [ ] **2.3** On failure: record test name, assertion, actual vs expected values, and whether the bug is in `R_refactored/` or in the test expectations.

---

## 3. Run Python unit tests

- [ ] **3.1** Install package:
  ```bash
  pip install -e ".[dev]"
  ```

- [ ] **3.2** Run all tests:
  ```bash
  pytest tests/ -v --tb=short
  ```

- [ ] **3.3** All 65 tests must pass (see `TESTS_PYTHON.md` for what each test checks).

- [ ] **3.4** Cross-validation only (needs reference CSVs from step 1):
  ```bash
  pytest tests/ -v -k "cross_validation or medium"
  ```

- [ ] **3.5** Unit tests only (no reference files needed):
  ```bash
  pytest tests/ -v -k "not cross_validation and not medium"
  ```

**SLURM shortcut:** `sbatch slurm/run_py_tests.slurm`

---

## 4. Investigate R vs Python discrepancies

When a cross-validation test fails, follow this procedure:

- [ ] **4.1** Identify the exact (signature, dataset, metric) triple that disagrees.
- [ ] **4.2** Check if it is a known acceptable difference (PCA sign, GaussianMixture tolerance). See `TESTS_PYTHON.md § test_compare_metrics.py` and `GUIDE.md § Key R/Python differences`.
- [ ] **4.3** Print intermediate values from both sides to find where the divergence starts:
  ```python
  result = compute_var(sigs, names_sigs, datasets, names_datasets)
  print(result["mean_sd_tables"]["compact_sig"]["dataset_A"])
  print(result["all_sd"]["compact_sig"]["dataset_A"].head(20))
  ```
- [ ] **4.4** Run the same step in R interactively:
  ```r
  source("R_refactored/eval_var_loc.R")
  r <- compute_var(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)
  r$radar_values[["compact_sig"]][["dataset_A"]]
  ```
- [ ] **4.5** Fix the Python bug in `pysigqc/`, re-run cross-validation, then propagate the same fix to `pysigqc_joblib/`.

Per-module pitfalls to check first:
- `eval_expr`: R's `na.omit()` on a matrix drops entire rows — Python must use `dropna(axis=0, how="any")`, not per-cell NA removal.
- `eval_compactness`: R's `cor()` handles constant columns per-pair — the `fill_diagonal` fix is needed in `spearman_matrix()`. See `OPTIMIZATIONS.md`.
- `eval_stan`: zero-variance genes must produce score 0, not NaN/Inf. See `AUDIT.md § BUG-2`.
- `compare_metrics`: PCA sign ambiguity — always compare absolute values for `rho_pca1_*` metrics.

---

## 5. Verify pysigqc vs pysigqc_joblib equivalence

- [ ] **5.1** Run equivalence tests:
  ```bash
  pytest tests/test_joblib_equivalence.py -v
  ```

- [ ] **5.2** All 11 tests must pass. Any failure is a bug in `pysigqc_joblib`. See `TESTS_PYTHON.md § test_joblib_equivalence.py` and `OPTIMIZATIONS.md` for which changes were made and where subtle differences can arise.

- [ ] **5.3** Verify parallel execution gives the same result as sequential:
  ```python
  from pysigqc_joblib.pipeline import run_pipeline
  r1 = run_pipeline(sigs, names_sigs, datasets, names_datasets, n_jobs=1)
  r2 = run_pipeline(sigs, names_sigs, datasets, names_datasets, n_jobs=-1)
  # radar tables should be identical
  ```

---

## 6. Performance benchmarks

### 6.1 Medium fixture — per-module timing

- [ ] For each module, compare `elapsed_seconds` between `pysigqc` and `pysigqc_joblib`:
  ```python
  import pysigqc.eval_compactness as ref
  import pysigqc_joblib.eval_compactness as opt
  r = ref.compute_compactness(sigs, names_sigs, datasets, names_datasets)
  o = opt.compute_compactness(sigs, names_sigs, datasets, names_datasets)
  print(f"speedup: {r['elapsed_seconds'] / o['elapsed_seconds']:.1f}x")
  ```

- [ ] Record results (expected: `eval_compactness` ~25x, others ~2–5x):

  | Module | pysigqc (s) | pysigqc_joblib (s) | Speedup |
  |--------|-------------|-------------------|---------|
  | eval_var | | | |
  | eval_expr | | | |
  | eval_compactness | | | |
  | eval_stan | | | |
  | compare_metrics | | | |

### 6.2 Scaling benchmark

- [ ] Generate synthetic datasets at increasing sizes and measure total pipeline time:

  | n_genes | n_samples | pysigqc (s) | pysigqc_joblib (s) | Speedup |
  |---------|-----------|-------------|-------------------|---------|
  | 100 | 50 | | | |
  | 1,000 | 1,000 | | | |
  | 10,000 | 10,000 | | | |
  | 20,000 | 100,000 | | | |

  ```python
  import numpy as np, pandas as pd
  def make_synthetic(n_genes, n_samples, seed=42):
      rng = np.random.default_rng(seed)
      arr = rng.normal(5, 2, size=(n_genes, n_samples))
      return pd.DataFrame(arr,
          index=[f"GENE{i:05d}" for i in range(n_genes)],
          columns=[f"S{i:06d}" for i in range(n_samples)])
  ```

### 6.3 Negative control parallel speedup

- [ ] Measure wall time vs `n_jobs`:

  | n_jobs | Time (s) | Speedup |
  |--------|----------|---------|
  | 1 | | 1.0x |
  | 4 | | |
  | 8 | | |
  | -1 | | |

---

## 7. BRCA paper data validation

- [ ] **7.1** Place BRCA data in `sigQC-master/input_files/` (expression CSVs + signature CSVs).
- [ ] **7.2** Run cross-validation job: `sbatch slurm/run_brca_cross_validation.slurm`
- [ ] **7.3** Inspect mismatches: `cat brca_cross_validation/py_output/brca_mismatches.csv`
- [ ] **7.4** For each mismatch: check if it is PCA sign, GaussianMixture tolerance, or a genuine discrepancy. Follow procedure in section 4 above.

---

## 8. Edge case stress tests

- [ ] **8.1** All-NA signature gene: set one sig gene to all NaN, verify no crash and NaN propagates correctly.
- [ ] **8.2** Single-gene signature: `compactness` should return `autocor_median=0` (no pairs), `compare_metrics` should skip PCA.
- [ ] **8.3** All-constant dataset: z-transform should return zeros, Spearman correlations NaN.
- [ ] **8.4** Signature with no genes in dataset: all modules should handle gracefully.
- [ ] **8.5** Large signature (500 genes): compactness should complete in reasonable time with the `corrcoef` optimisation.

---

## 9. Reproducibility

- [ ] **9.1** `pysigqc` is deterministic (no randomness): verify two runs give identical `radar_values`.
- [ ] **9.2** `pysigqc_joblib` negative control with `n_jobs=1` and same seed: identical `metrics_table`.
- [ ] **9.3** `pysigqc_joblib` negative control with `n_jobs=-1` and same seed: identical results across two parallel runs (uses `SeedSequence.spawn()`).
- [ ] **9.4** Cross-platform: run same fixture on macOS (local) and Linux (cluster), verify values match within `rtol=1e-12`.

---

## 10. Cleanup

- [ ] **10.1** After benchmarks, fill in the blank tables in sections 6.1–6.3 and update `OPTIMIZATIONS.md`.
- [ ] **10.2** Update `GUIDE.md` with final test count and any new edge cases discovered.
- [ ] **10.3** Verify `git status` is clean (no untracked source files).
- [ ] **10.4** Final full test run: `pytest tests/ -v --tb=short` — all 65 tests pass.
