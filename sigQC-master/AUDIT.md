# sigQC Code Audit Report

Audit performed on the sigQC R package (v0.1.24) as part of the code review and
Python porting project. All issues are annotated in the source files with `# BUG:` or
`# ISSUE:` tags.

---

## Confirmed Bugs

### BUG-1: Discrete gene sampling uses `runif()` instead of `sample()` [CRITICAL]
- **File:** `R/sigsQcNegativeControl.R`, line 176
- **Code:** `random.index.vector = stats::runif(min=1, max=data.matrix.nrows, n=len)`
- **Problem:** `runif()` returns continuous values. When used as row indices, R truncates to integer, causing biased sampling: gene at index 1 has near-zero probability, last gene gets extra probability, and genes can be selected more than once (no uniqueness).
- **Fix:** `sample(1:data.matrix.nrows, size=len, replace=FALSE)`

### BUG-2: Division by zero in z-transform for zero-variance genes [CRITICAL]
- **Files:** `R/eval_stan_loc.R` line 54, `R/eval_struct_loc.R` lines 208/265/345/399
- **Code:** `(expr - mean(expr)) / sd(expr)` where `sd()` can return 0
- **Problem:** Genes with constant expression across all samples have `sd() == 0`. Division produces `NaN`/`Inf` that propagates silently through `na.omit()`, changing the effective signature size without warning.
- **Fix:** Guard: `if (sd(...) == 0) { z_score <- 0 }`

### BUG-3: Index bug when removing short signatures [CRITICAL]
- **File:** `R/make_all_plots.R`, line 162
- **Code:** `gene_sigs_list[names_sigs[k]] <- NULL`
- **Problem:** Uses loop variable `k` (which holds the last value from the preceding loop = `length(names_sigs)`) instead of iterating over `sigs_to_remove_ind`. Only the last signature is removed, not all flagged ones.
- **Fix:** `gene_sigs_list[names_sigs[sigs_to_remove_ind]] <- NULL`

### BUG-4: Mclust model name labels E/V swapped [MEDIUM]
- **File:** `R/compare_metrics_loc.R`, lines 466-469 (and repeated at lines 496-499, 525-528)
- **Problem:** When `modelName == 'V'`, the output says "equal variance (E)" and vice versa. The correct mapping is: E = Equal variance, V = Variable variance.
- **Fix:** Swap the condition labels.

### BUG-5: Mean mixture model fitted on median scores (copy-paste error) [MEDIUM]
- **File:** `R/compare_metrics_loc.R`, line 486
- **Code:** `mclust::Mclust(med_scores, G=1:max_clusters)` — stored as `[['mean']]`
- **Problem:** The "mean" mixture model is fitted on `med_scores` instead of `mean_scores`, making it a duplicate of the median model.
- **Fix:** Replace `med_scores` with `mean_scores`.

### BUG-6: `na.omit(t(...))` removes samples instead of genes [MEDIUM]
- **File:** `R/compare_metrics_loc.R`, line 56 (and lines 223, 437, 592)
- **Code:** `stats::prcomp(stats::na.omit(t(data.matrix[inter,])), retx=T)`
- **Problem:** `t()` transposes the matrix so rows = samples, columns = genes. `na.omit()` then removes entire samples that have any NA in their gene values, rather than removing genes with NAs.
- **Impact:** With even one missing gene value, entire samples are dropped from PCA.

### BUG-7: Operator precedence in length check [MEDIUM]
- **File:** `R/eval_struct_loc.R`, line 61
- **Code:** `if(length(rows_needed >0))`
- **Problem:** Evaluates `rows_needed > 0` first (comparing character to numeric), then takes `length()`. Works by accident but for the wrong reason.
- **Fix:** `if(length(rows_needed) > 0)`

### BUG-8: Missing gene padding with min expression value [MEDIUM]
- **Files:** `R/eval_struct_loc.R` line 62, `R/eval_compactness_loc.R` line 62
- **Code:** `matrix(min(sig_scores), nrow=length(rows_needed), ...)`
- **Problem:** Padding missing genes with the dataset minimum creates rows with identical values across all samples, introducing artificial zero-variance rows that distort clustering and autocorrelation.
- **Fix:** Pad with `NA` and use `na.rm=TRUE` downstream.

---

## Theoretical / Methodological Issues

### ISSUE-1: Binarization threshold after z-transform is questionable
- **File:** `R/eval_struct_loc.R`, lines 212/269/348/402
- **Threshold:** `min + (max - min) / 2` (midpoint of the range)
- **Problem:** For z-scored data, this depends on outliers (min and max) rather than any statistical property. A more principled approach: binarize at 0 (the mean of z-scores) or use Otsu's method.

### ISSUE-2: `abs()` applied to radar chart metrics loses directionality
- **File:** `R/make_radar_chart_loc.R`, line 44
- **Problem:** `abs(radar_plot_mat)` converts all negative correlations to positive. Negative autocorrelation (anti-correlation between signature genes) is a meaningful QC signal, but is indistinguishable from positive autocorrelation on the chart.

### ISSUE-3: Radar chart area formula is incorrect
- **File:** `R/make_radar_chart_loc.R`, line 106
- **Formula used:** `sum(r_i * r_{i+1}) / n`
- **Correct formula:** `0.5 * sum(r_i * r_{i+1} * sin(2*pi/n))`
- **Impact:** Missing `sin()` factor means computed values are a dimensionless proxy, not actual areas. Ordinal comparisons (larger = better) still hold.

### ISSUE-4: No multiple testing correction
- **File:** `R/compare_metrics_loc.R` (throughout)
- **Problem:** Many Spearman correlations are computed (mean-median, mean-PCA1, PCA1-median, GSVA-ssGSEA, ssGSEA-PLAGE, PLAGE-GSVA per signature per dataset) but raw rho values are reported without p-values or FDR correction.

### ISSUE-5: GSVA/ssGSEA/PLAGE treated interchangeably
- **File:** `R/compare_metrics_loc.R`, Section 2
- **Problem:** These methods have fundamentally different statistical assumptions (nonparametric kernel density vs. rank-based vs. SVD). Pairwise correlations assess concordance but no guidance is given on which is more appropriate.

### ISSUE-6: BCCC biclustering parameters are arbitrary
- **File:** `R/eval_struct_loc.R`, lines 231/291/366/420
- **Parameters:** `delta=1, alpha=1.5, number=50`
- **Problem:** These appear to be BCCC defaults with no domain justification. The `num_cols_chosen` variable computed nearby is dead code left from a removed `BCXmotifs` call.

### ISSUE-7: Permutation strategy destroys sample-level structure
- **File:** `R/sigsQcNegativeControl.R`, line 285
- **Design:** Each sample column gets an independent random permutation of gene labels.
- **Question:** This destroys both inter-gene correlation AND sample-level patterns. If the goal is only to test inter-gene structure, a single permutation applied to all samples would preserve sample correlations while still breaking gene identity.

---

## Code Quality Issues

### DRY Violation: ~40% code duplication
- `R/compute_without_plots.R` (398 lines) duplicates all `eval_*_loc` functions without plotting code.
- `R/eval_struct_loc.R` repeats the z-transform → binarize → BCCC biclustering pipeline 4 times.
- `R/compare_metrics_loc.R` recomputes mean/median/PCA1 scores in every section.

### Broad `tryCatch` blocks mask errors
- `R/make_all_plots.R` wraps each pipeline step in `tryCatch` that logs but doesn't halt, so failures in early steps (e.g., eval_var_loc) are silently logged and subsequent steps may produce incomplete results.

### Magic numbers throughout
- PDF sizing: 25/50/75 element thresholds in boxplot.matrix2
- Font scaling: `5/log10(h)`, `3*10/max_title_length`, `4*10/max_title_length`
- BCCC parameters: `delta=1, alpha=1.5, number=50`
- Quantile thresholds: 0.025, 0.25, 0.5, 0.75, 0.975

### Parameter naming inconsistency
- `sigsQcNegativeControl.R`: outer function uses `numResampling`, inner uses `rePower` — same concept, different names.
