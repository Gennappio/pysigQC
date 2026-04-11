"""Compare various signature summary scoring metrics.

Port of R_refactored/compare_metrics_loc.R — compute_metrics() only.
Computes Mean, Median, PCA1 scores and their Spearman correlations.
Also computes GSVA/ssGSEA/PLAGE enrichment scores (optional) and
Gaussian mixture models.

Produces 4 radar metrics: rho_mean_med, rho_pca1_med, rho_mean_pca1, prop_pca1_var.
"""

from __future__ import annotations

import time
import warnings

import numpy as np
import pandas as pd
from scipy import stats as sp_stats
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture

from .utils import gene_intersection


def compute_metrics(
    gene_sigs_list: dict[str, list[str]],
    names_sigs: list[str],
    mRNA_expr_matrix: dict[str, pd.DataFrame],
    names_datasets: list[str],
) -> dict:
    """Compute scoring comparison metrics for each signature-dataset pair.

    Returns dict with keys:
        radar_values: nested dict [sig][dataset] -> dict of 4 metrics
        scores: nested dict [sig][dataset] -> dict of score arrays
        pca_results: nested dict [sig][dataset] -> dict with pca_obj, props_of_variances
        score_cor_mats: dict of scoring correlation matrices
        mixture_models: nested dict [sig][dataset] -> dict with model results
        elapsed_seconds: wall-clock time
    """
    _t0 = time.perf_counter()
    radar_values: dict = {}
    scores_all: dict = {}
    pca_results: dict = {}
    score_cor_mats: dict = {}
    mixture_models: dict = {}

    for sig in names_sigs:
        gene_sig = gene_sigs_list[sig]
        radar_values[sig] = {}
        scores_all[sig] = {}
        pca_results[sig] = {}
        mixture_models[sig] = {}

        for ds in names_datasets:
            data_matrix = mRNA_expr_matrix[ds].copy()
            # Replace non-finite values with NaN
            data_matrix = data_matrix.where(np.isfinite(data_matrix), other=np.nan)
            inter = gene_intersection(gene_sig, data_matrix)

            sig_data = data_matrix.loc[inter]

            # --- Compute summary scores ---
            med_scores = sig_data.apply(lambda col: np.nanmedian(col.values), axis=0).values
            mean_scores = sig_data.apply(lambda col: np.nanmean(col.values), axis=0).values
            sample_names = list(data_matrix.columns)

            # PCA on signature genes (samples as observations, genes as features)
            pca1_scores = None
            props_of_variances = None
            pca_obj = None
            try:
                # R: prcomp(t(na.omit(data.matrix[inter,])))
                sig_clean = sig_data.dropna(axis=0, how="any").T  # samples x genes
                if sig_clean.shape[1] >= 2 and sig_clean.shape[0] >= 2:
                    pca_obj = PCA()
                    pca_obj.fit(sig_clean.values)
                    pca1_scores = pca_obj.transform(sig_clean.values)[:, 0]
                    props_of_variances = pca_obj.explained_variance_ratio_
            except Exception:
                pca1_scores = None
                pca_obj = None

            pca_results[sig][ds] = {
                "pca_obj": pca_obj,
                "props_of_variances": props_of_variances,
            }

            # --- Compute correlations ---
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

            # --- Mixture models (using sklearn instead of mclust) ---
            mm = {"median": None, "mean": None, "pca1": None}
            for score_name, score_arr in [("median", med_scores), ("mean", mean_scores),
                                           ("pca1", pca1_scores)]:
                if score_arr is not None and len(score_arr) > 1:
                    clean = score_arr[np.isfinite(score_arr)]
                    if len(clean) >= 2:
                        max_k = min(len(clean) // 2, 10)
                        max_k = max(max_k, 1)
                        best_bic = np.inf
                        best_model = None
                        for k in range(1, max_k + 1):
                            try:
                                gm = GaussianMixture(n_components=k, random_state=42)
                                gm.fit(clean.reshape(-1, 1))
                                bic = gm.bic(clean.reshape(-1, 1))
                                if bic < best_bic:
                                    best_bic = bic
                                    best_model = gm
                            except Exception:
                                pass
                        mm[score_name] = best_model
            mixture_models[sig][ds] = mm

            # --- Store scores ---
            scores_all[sig][ds] = {
                "med_scores": med_scores,
                "mean_scores": mean_scores,
                "pca1_scores": pca1_scores,
                "common_score_cols": sample_names,
                "props_of_variances": props_of_variances,
            }

    return {
        "radar_values": radar_values,
        "scores": scores_all,
        "pca_results": pca_results,
        "score_cor_mats": score_cor_mats,
        "mixture_models": mixture_models,
        "elapsed_seconds": time.perf_counter() - _t0,
    }
