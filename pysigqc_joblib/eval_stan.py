"""Evaluate standardization effects — vectorized numpy version.

Same API and return types as pysigqc.eval_stan.compute_stan().
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd
from scipy import stats as sp_stats

from ._core import to_numpy, gene_indices, z_transform_matrix, nanmedian_cols
from .utils import gene_intersection


def compute_stan(
    gene_sigs_list: dict[str, list[str]],
    names_sigs: list[str],
    mRNA_expr_matrix: dict[str, pd.DataFrame],
    names_datasets: list[str],
) -> dict:
    """Compute standardization comparison metrics."""
    _t0 = time.perf_counter()
    radar_values: dict = {}
    med_scores_all: dict = {}
    z_transf_scores_all: dict = {}

    # Pre-convert datasets to numpy once
    ds_cache: dict = {}
    for ds in names_datasets:
        ds_cache[ds] = to_numpy(mRNA_expr_matrix[ds])

    for sig in names_sigs:
        gene_sig = gene_sigs_list[sig]
        radar_values[sig] = {}
        med_scores_all[sig] = {}
        z_transf_scores_all[sig] = {}

        for ds in names_datasets:
            arr, row_names, col_names = ds_cache[ds]
            inter = gene_intersection(gene_sig, mRNA_expr_matrix[ds])
            sig_idx = gene_indices(inter, row_names)

            sig_arr = arr[sig_idx].copy()  # (n_genes, n_samples)

            # Vectorized z-transform all genes at once
            z_arr = z_transform_matrix(sig_arr)

            # Column-wise median (across genes for each sample)
            med_scores = nanmedian_cols(sig_arr)
            z_scores = nanmedian_cols(z_arr)

            # Spearman correlation between raw and z-transformed
            rho, _ = sp_stats.spearmanr(med_scores, z_scores)

            radar_values[sig][ds] = {"standardization_comp": float(rho)}
            med_scores_all[sig][ds] = med_scores
            z_transf_scores_all[sig][ds] = z_scores

    return {
        "radar_values": radar_values,
        "med_scores": med_scores_all,
        "z_transf_scores": z_transf_scores_all,
        "elapsed_seconds": time.perf_counter() - _t0,
    }
