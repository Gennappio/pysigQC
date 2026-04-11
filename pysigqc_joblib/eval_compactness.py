"""Evaluate compactness via autocorrelation — vectorized with fast Spearman.

Same API and return types as pysigqc.eval_compactness.compute_compactness().

Key optimization: instead of O(n_genes^2) individual spearmanr() calls,
we rank all genes once and use np.corrcoef on the rank matrix (single BLAS call).
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd

from ._core import to_numpy, gene_indices, spearman_matrix
from .utils import gene_intersection


def compute_compactness(
    gene_sigs_list: dict[str, list[str]],
    names_sigs: list[str],
    mRNA_expr_matrix: dict[str, pd.DataFrame],
    names_datasets: list[str],
) -> dict:
    """Compute compactness metrics for each signature-dataset pair."""
    _t0 = time.perf_counter()
    radar_values: dict = {}
    autocor_matrices: dict = {}

    # Pre-convert datasets to numpy once
    ds_cache: dict = {}
    for ds in names_datasets:
        ds_cache[ds] = to_numpy(mRNA_expr_matrix[ds])

    for sig in names_sigs:
        gene_sig = gene_sigs_list[sig]
        radar_values[sig] = {}
        autocor_matrices[sig] = {}

        for ds in names_datasets:
            arr, row_names, col_names = ds_cache[ds]
            inter = gene_intersection(gene_sig, mRNA_expr_matrix[ds])

            # Get signature gene rows, drop genes with any NaN (matches R's na.omit)
            sig_idx = gene_indices(inter, row_names)
            sig_arr = arr[sig_idx]  # (n_sig_genes, n_samples)

            # Drop rows with any NaN
            clean_mask = ~np.isnan(sig_arr).any(axis=1)
            sig_clean = sig_arr[clean_mask]
            genes_present = [inter[i] for i in range(len(inter)) if clean_mask[i]]

            # Spearman correlation matrix — fast path uses rank + corrcoef
            autocors = spearman_matrix(sig_clean)

            autocor_df = pd.DataFrame(autocors, index=genes_present, columns=genes_present)
            autocor_matrices[sig][ds] = autocor_df

            if autocors.shape[0] > 1:
                autocor_median = float(np.nanmedian(autocors))
            else:
                autocor_median = 0.0

            radar_values[sig][ds] = {"autocor_median": autocor_median}

    return {
        "radar_values": radar_values,
        "autocor_matrices": autocor_matrices,
        "rank_product_tables": {},
        "elapsed_seconds": time.perf_counter() - _t0,
    }
