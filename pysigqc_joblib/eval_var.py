"""Evaluate expression variability — vectorized numpy version.

Same API and return types as pysigqc.eval_var.compute_var().
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd

from ._core import to_numpy, gene_indices, nanstd_rows, nanmean_rows, skew_rows
from .utils import gene_intersection


def compute_var(
    gene_sigs_list: dict[str, list[str]],
    names_sigs: list[str],
    mRNA_expr_matrix: dict[str, pd.DataFrame],
    names_datasets: list[str],
) -> dict:
    """Compute variability metrics for each signature-dataset pair."""
    _t0 = time.perf_counter()
    radar_values: dict = {}
    mean_sd_tables: dict = {}
    all_sd: dict = {}
    all_mean: dict = {}
    inter_genes: dict = {}

    # Pre-compute per-dataset numpy arrays and stats (done once per dataset)
    ds_cache: dict = {}
    for ds in names_datasets:
        arr, row_names, col_names = to_numpy(mRNA_expr_matrix[ds])
        sd_all = nanstd_rows(arr, ddof=1)
        mean_all = nanmean_rows(arr)
        cv_all = sd_all / mean_all  # may have inf/nan for mean=0
        ds_cache[ds] = {
            "arr": arr, "row_names": row_names, "col_names": col_names,
            "sd_all": sd_all, "mean_all": mean_all, "cv_all": cv_all,
        }

    for sig in names_sigs:
        gene_sig = gene_sigs_list[sig]
        radar_values[sig] = {}
        mean_sd_tables[sig] = {}
        all_sd[sig] = {}
        all_mean[sig] = {}
        inter_genes[sig] = {}

        for ds in names_datasets:
            cache = ds_cache[ds]
            arr = cache["arr"]
            row_names = cache["row_names"]
            sd_all = cache["sd_all"]
            mean_all = cache["mean_all"]
            cv_all = cache["cv_all"]

            inter = gene_intersection(gene_sig, mRNA_expr_matrix[ds])
            inter_genes[sig][ds] = inter
            sig_idx = gene_indices(inter, row_names)

            all_sd[sig][ds] = pd.Series(sd_all, index=row_names)
            all_mean[sig][ds] = pd.Series(mean_all, index=row_names)

            # Signature subset
            sd_sig = sd_all[sig_idx]
            sd_sig_clean = sd_sig[~np.isnan(sd_sig)]
            sd_all_clean = sd_all[~np.isnan(sd_all)]

            # Metric 1: sd_median_ratio
            med_sig_sd = np.median(sd_sig_clean) if len(sd_sig_clean) > 0 else 0.0
            med_all_sd = np.median(sd_all_clean) if len(sd_all_clean) > 0 else 0.0
            denom = med_all_sd + med_sig_sd
            sd_median_ratio = med_sig_sd / denom if denom != 0 else 0.0

            # Metric 2: abs_skewness_ratio
            mean_sig = mean_all[sig_idx]
            mean_sig_clean = mean_sig[~np.isnan(mean_sig)]
            mean_all_clean = mean_all[~np.isnan(mean_all)]

            skew_sig = abs(float(skew_rows(mean_sig_clean.reshape(1, -1))[0])) if len(mean_sig_clean) > 2 else 0.0
            skew_all = abs(float(skew_rows(mean_all_clean.reshape(1, -1))[0])) if len(mean_all_clean) > 2 else 0.0
            denom = skew_all + skew_sig
            abs_skewness_ratio = skew_sig / denom if denom != 0 else 0.0

            # Mean/SD table
            mean_sd_tables[sig][ds] = pd.DataFrame(
                {"Mean": mean_all[sig_idx], "SD": sd_all[sig_idx]},
                index=inter,
            )

            # CV for signature genes
            cv_sig = cv_all[sig_idx]
            cv_sig_clean = cv_sig[np.isfinite(cv_sig)]
            cv_all_clean = cv_all[np.isfinite(cv_all)]

            # Metrics 3-5: proportion in top quantiles
            if len(cv_all_clean) > 0 and len(cv_sig_clean) > 0:
                q90, q75, q50 = np.nanquantile(cv_all_clean, [0.9, 0.75, 0.5])
                n_sig = len(cv_sig_clean)
                prop_top_10 = float(np.sum(cv_sig_clean >= q90)) / n_sig
                prop_top_25 = float(np.sum(cv_sig_clean >= q75)) / n_sig
                prop_top_50 = float(np.sum(cv_sig_clean >= q50)) / n_sig
            else:
                prop_top_10 = prop_top_25 = prop_top_50 = 0.0

            # Metric 6: coeff_of_var_ratio
            med_cv_sig = abs(np.nanmedian(cv_sig_clean)) if len(cv_sig_clean) > 0 else 0.0
            med_cv_all = abs(np.nanmedian(cv_all_clean)) if len(cv_all_clean) > 0 else 0.0
            denom = med_cv_all + med_cv_sig
            coeff_of_var_ratio = med_cv_sig / denom if denom != 0 else 0.0

            radar_values[sig][ds] = {
                "sd_median_ratio": sd_median_ratio,
                "abs_skewness_ratio": abs_skewness_ratio,
                "prop_top_10_percent": prop_top_10,
                "prop_top_25_percent": prop_top_25,
                "prop_top_50_percent": prop_top_50,
                "coeff_of_var_ratio": coeff_of_var_ratio,
            }

    return {
        "radar_values": radar_values,
        "mean_sd_tables": mean_sd_tables,
        "all_sd": all_sd,
        "all_mean": all_mean,
        "inter": inter_genes,
        "elapsed_seconds": time.perf_counter() - _t0,
    }
