"""Evaluate expression level properties — vectorized numpy version.

Same API and return types as pysigqc.eval_expr.compute_expr().
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd

from ._core import to_numpy, gene_indices, rows_without_nan
from .utils import gene_intersection


def compute_expr(
    gene_sigs_list: dict[str, list[str]],
    names_sigs: list[str],
    mRNA_expr_matrix: dict[str, pd.DataFrame],
    names_datasets: list[str],
    thresholds: dict[str, float] | None = None,
) -> dict:
    """Compute expression-level metrics for each signature-dataset pair."""
    _t0 = time.perf_counter()
    radar_values: dict = {}
    na_proportions: dict = {}
    expr_proportions: dict = {}

    # Pre-compute per-dataset numpy arrays and thresholds
    ds_cache: dict = {}
    computed_thresholds: dict[str, float] = {}

    for ds in names_datasets:
        arr, row_names, col_names = to_numpy(mRNA_expr_matrix[ds])
        ds_cache[ds] = {"arr": arr, "row_names": row_names, "col_names": col_names}

        # Threshold: median of entire matrix after dropping rows with any NaN
        if thresholds is None:
            clean_mask = rows_without_nan(arr)
            clean_vals = arr[clean_mask].ravel()
            computed_thresholds[ds] = float(np.median(clean_vals)) if len(clean_vals) > 0 else 0.0
        elif isinstance(thresholds, dict):
            computed_thresholds[ds] = thresholds[ds]

    if thresholds is not None and not isinstance(thresholds, dict):
        computed_thresholds = dict(zip(names_datasets, thresholds))

    # --- NA proportion + expression proportion ---
    for sig in names_sigs:
        gene_sig = gene_sigs_list[sig]
        radar_values[sig] = {}
        na_proportions[sig] = {}
        expr_proportions[sig] = {}

        for ds in names_datasets:
            cache = ds_cache[ds]
            arr = cache["arr"]
            row_names = cache["row_names"]
            n_samples = arr.shape[1]

            inter = gene_intersection(gene_sig, mRNA_expr_matrix[ds])
            sig_idx = gene_indices(inter, row_names)

            # Vectorized NA proportions per gene
            sig_arr = arr[sig_idx]  # (n_sig_genes, n_samples)
            na_counts = np.isnan(sig_arr).sum(axis=1)
            gene_na_props_vals = na_counts / n_samples

            # Missing genes get NA proportion = 1.0
            missing_genes = [g for g in gene_sig if g not in set(row_names)]
            all_gene_names = list(inter) + missing_genes
            all_na_props = np.concatenate([gene_na_props_vals, np.ones(len(missing_genes))])

            # Sort descending
            sort_idx = np.argsort(-all_na_props)
            na_proportions[sig][ds] = pd.Series(
                all_na_props[sort_idx],
                index=[all_gene_names[i] for i in sort_idx],
            )

            # Radar metric: median proportion of non-NA
            radar_values[sig][ds] = {
                "med_prop_na": float(np.median(1 - all_na_props)),
            }

            # --- Expression proportions ---
            thresh = computed_thresholds[ds]

            # Vectorized: proportion of samples above threshold per gene
            # NaN comparisons → False, so we need to handle NA propagation
            below = sig_arr < thresh  # NaN → False
            has_na = np.isnan(sig_arr).any(axis=1)
            props = 1.0 - below.sum(axis=1) / n_samples
            # Match R: genes with any NA get NaN proportion
            props[has_na] = np.nan

            # Missing genes get 0.0
            all_props = np.concatenate([props, np.zeros(len(missing_genes))])

            # Sort ascending
            sort_idx = np.argsort(all_props)  # NaN goes to end
            expr_proportions[sig][ds] = pd.Series(
                all_props[sort_idx],
                index=[all_gene_names[i] for i in sort_idx],
            )

            # R's sort() removes NAs, then median
            radar_values[sig][ds]["med_prop_above_med"] = float(np.nanmedian(all_props))

    return {
        "radar_values": radar_values,
        "na_proportions": na_proportions,
        "expr_proportions": expr_proportions,
        "thresholds": computed_thresholds,
        "elapsed_seconds": time.perf_counter() - _t0,
    }
