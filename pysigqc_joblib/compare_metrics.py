"""Compare scoring metrics — IncrementalPCA for large datasets.

Same API and return types as pysigqc.compare_metrics.compute_metrics().

Key optimization: uses IncrementalPCA for datasets with >50k samples to avoid
loading entire (n_samples x n_genes) matrix into memory at once.
"""

from __future__ import annotations

import time
import warnings

import numpy as np
import pandas as pd
from scipy import stats as sp_stats
from sklearn.decomposition import PCA, IncrementalPCA
from sklearn.mixture import GaussianMixture

from ._core import to_numpy, gene_indices, rows_without_nan, nanmedian_cols, nanmean_rows
from .utils import gene_intersection

# Import enrichment functions from reference implementation
from pysigqc.compare_metrics import _gsva_score, _ssgsea_score, _plage_score

INCREMENTAL_THRESHOLD = 50_000


def _fit_pca(sig_clean_T: np.ndarray, n_genes: int):
    """Fit PCA, using IncrementalPCA for large sample counts.

    sig_clean_T: (n_samples, n_genes) — already transposed and NaN-free.
    Returns (pca1_scores, explained_variance_ratio, pca_obj).
    """
    n_samples = sig_clean_T.shape[0]

    if n_samples > INCREMENTAL_THRESHOLD:
        n_components = min(n_genes, 10)
        ipca = IncrementalPCA(n_components=n_components)
        chunk_size = max(n_components + 1, 10_000)
        for start in range(0, n_samples, chunk_size):
            end = min(start + chunk_size, n_samples)
            if (end - start) > n_components:
                ipca.partial_fit(sig_clean_T[start:end])
        pca1_scores = ipca.transform(sig_clean_T)[:, 0]
        return pca1_scores, ipca.explained_variance_ratio_, ipca
    else:
        pca_obj = PCA()
        pca_obj.fit(sig_clean_T)
        pca1_scores = pca_obj.transform(sig_clean_T)[:, 0]
        return pca1_scores, pca_obj.explained_variance_ratio_, pca_obj


def compute_metrics(
    gene_sigs_list: dict[str, list[str]],
    names_sigs: list[str],
    mRNA_expr_matrix: dict[str, pd.DataFrame],
    names_datasets: list[str],
    compute_enrichment: bool = True,
    fit_mixture: bool = True,
) -> dict:
    """Compute scoring comparison metrics for each signature-dataset pair.

    When ``compute_enrichment=False`` the GSVA/ssGSEA/PLAGE scores are skipped;
    when ``fit_mixture=False`` the Gaussian mixture models are skipped. Both
    outputs are only consumed by plot_metrics and by the radar-chart pipeline
    indirectly, so the negative/permutation control callers can opt out.
    """
    _t0 = time.perf_counter()
    radar_values: dict = {}
    scores_all: dict = {}
    pca_results: dict = {}
    score_cor_mats: dict = {}
    mixture_models: dict = {}
    enrichment_scores: dict = {}

    # Pre-convert datasets
    ds_cache: dict = {}
    for ds in names_datasets:
        arr, row_names, col_names = to_numpy(mRNA_expr_matrix[ds])
        # Replace non-finite with NaN
        arr = np.where(np.isfinite(arr), arr, np.nan)
        ds_cache[ds] = (arr, row_names, col_names)

    for sig in names_sigs:
        gene_sig = gene_sigs_list[sig]
        radar_values[sig] = {}
        scores_all[sig] = {}
        pca_results[sig] = {}
        mixture_models[sig] = {}
        enrichment_scores[sig] = {}

        for ds in names_datasets:
            arr, row_names, col_names = ds_cache[ds]
            inter = gene_intersection(gene_sig, mRNA_expr_matrix[ds])
            sig_idx = gene_indices(inter, row_names)

            sig_arr = arr[sig_idx]  # (n_genes, n_samples)

            # Vectorized median/mean scores per sample
            with np.errstate(invalid="ignore"):
                med_scores = np.nanmedian(sig_arr, axis=0)
                mean_scores = np.nanmean(sig_arr, axis=0)

            # PCA on clean (no-NaN) gene rows
            pca1_scores = None
            props_of_variances = None
            pca_obj = None
            try:
                clean_mask = rows_without_nan(sig_arr)
                sig_clean = sig_arr[clean_mask]  # (n_clean_genes, n_samples)
                if sig_clean.shape[0] >= 2 and sig_clean.shape[1] >= 2:
                    sig_clean_T = sig_clean.T  # (n_samples, n_genes)
                    pca1_scores, props_of_variances, pca_obj = _fit_pca(
                        sig_clean_T, sig_clean.shape[0]
                    )
            except (np.linalg.LinAlgError, ValueError) as e:
                # Narrow to numerical failures R's prcomp equivalent can raise.
                # Matches R_refactored/compare_metrics_loc.R:52.
                warnings.warn(
                    f"PCA failed for sig={sig!r} ds={ds!r}: {type(e).__name__}: {e}",
                    RuntimeWarning, stacklevel=2,
                )
                pca1_scores = None
                pca_obj = None

            pca_results[sig][ds] = {
                "pca_obj": pca_obj,
                "props_of_variances": props_of_variances,
            }

            # --- Compute enrichment scores (GSVA, ssGSEA, PLAGE) ---
            gsva_scores = ssgsea_scores = plage_scores = None
            if compute_enrichment and len(inter) >= 2:
                sig_idx_list = list(sig_idx)
                if len(sig_idx_list) >= 2:
                    gsva_scores = _gsva_score(arr, sig_idx_list)
                    ssgsea_scores = _ssgsea_score(arr, sig_idx_list)
                    plage_scores = _plage_score(arr, sig_idx_list)

            enrichment_scores[sig][ds] = {
                "gsva": gsva_scores,
                "ssgsea": ssgsea_scores,
                "plage": plage_scores,
            }

            # Spearman correlations
            rho_mean_med = 0.0
            rho_mean_pca1 = 0.0
            rho_pca1_med = 0.0

            if len(med_scores) > 1 and len(mean_scores) > 1:
                rho_mean_med, _ = sp_stats.spearmanr(med_scores, mean_scores)

            if pca1_scores is not None and len(pca1_scores) > 1:
                rho_mean_pca1, _ = sp_stats.spearmanr(mean_scores, pca1_scores)
                rho_pca1_med, _ = sp_stats.spearmanr(pca1_scores, med_scores)

            prop_pca1_var = 0.0
            if props_of_variances is not None and len(props_of_variances) > 0:
                prop_pca1_var = float(props_of_variances[0])

            radar_values[sig][ds] = {
                "rho_mean_med": float(rho_mean_med),
                "rho_pca1_med": float(rho_pca1_med),
                "rho_mean_pca1": float(rho_mean_pca1),
                "prop_pca1_var": prop_pca1_var,
            }

            # Mixture models
            mm = {"median": None, "mean": None, "pca1": None}
            if fit_mixture:
                for score_name, score_arr in [("median", med_scores), ("mean", mean_scores),
                                               ("pca1", pca1_scores)]:
                    if score_arr is not None and len(score_arr) > 1:
                        clean = score_arr[np.isfinite(score_arr)]
                        if len(clean) >= 2:
                            max_k = min(len(clean) // 2, 10)
                            max_k = max(max_k, 1)
                            best_bic = np.inf
                            best_model = None
                            bic_values = []
                            for k in range(1, max_k + 1):
                                try:
                                    gm = GaussianMixture(n_components=k, random_state=42)
                                    gm.fit(clean.reshape(-1, 1))
                                    bic = gm.bic(clean.reshape(-1, 1))
                                    bic_values.append((k, bic))
                                    if bic < best_bic:
                                        best_bic = bic
                                        best_model = gm
                                except (ValueError, np.linalg.LinAlgError) as e:
                                    # BIC sweep: individual k-value fits can fail on
                                    # degenerate data; other k values may still succeed.
                                    warnings.warn(
                                        f"GaussianMixture(n_components={k}) failed for "
                                        f"sig={sig!r} ds={ds!r} score={score_name!r}: "
                                        f"{type(e).__name__}: {e}",
                                        RuntimeWarning, stacklevel=2,
                                    )
                            mm[score_name] = {
                                "best_model": best_model,
                                "bic_values": bic_values,
                                "best_k": best_model.n_components if best_model else None,
                            }
            mixture_models[sig][ds] = mm

            # --- Build scoring correlation matrix ---
            score_cols = {"Mean": med_scores, "Median": mean_scores}
            if pca1_scores is not None:
                score_cols["PCA1"] = pca1_scores
            # Add enrichment scores if available
            if gsva_scores is not None and np.isfinite(gsva_scores).any():
                score_cols["GSVA"] = gsva_scores
            if ssgsea_scores is not None and np.isfinite(ssgsea_scores).any():
                score_cols["ssGSEA"] = ssgsea_scores
            if plage_scores is not None and np.isfinite(plage_scores).any():
                score_cols["PLAGE"] = plage_scores

            if len(score_cols) >= 2:
                score_df = pd.DataFrame(score_cols)
                cor_mat = score_df.corr(method="spearman")
                score_cor_mats[f"{ds}_{sig}"] = cor_mat

            scores_all[sig][ds] = {
                "med_scores": med_scores,
                "mean_scores": mean_scores,
                "pca1_scores": pca1_scores,
                "gsva_scores": gsva_scores,
                "ssgsea_scores": ssgsea_scores,
                "plage_scores": plage_scores,
                "common_score_cols": col_names,
                "props_of_variances": props_of_variances,
            }

    return {
        "radar_values": radar_values,
        "scores": scores_all,
        "pca_results": pca_results,
        "score_cor_mats": score_cor_mats,
        "mixture_models": mixture_models,
        "enrichment_scores": enrichment_scores,
        "elapsed_seconds": time.perf_counter() - _t0,
    }
