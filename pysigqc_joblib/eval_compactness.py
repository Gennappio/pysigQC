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

# Import RankProd from reference implementation
from pysigqc.eval_compactness import _compute_rank_product


def compute_compactness(
    gene_sigs_list: dict[str, list[str]],
    names_sigs: list[str],
    mRNA_expr_matrix: dict[str, pd.DataFrame],
    names_datasets: list[str],
    compute_rank_product: bool = True,
    n_permutations: int = 100,
) -> dict:
    """Compute compactness metrics for each signature-dataset pair."""
    _t0 = time.perf_counter()
    radar_values: dict = {}
    autocor_matrices: dict = {}
    rank_product_tables: dict = {}
    gene_median_autocor: dict = {}

    # Pre-convert datasets to numpy once
    ds_cache: dict = {}
    for ds in names_datasets:
        ds_cache[ds] = to_numpy(mRNA_expr_matrix[ds])

    for sig in names_sigs:
        gene_sig = gene_sigs_list[sig]
        radar_values[sig] = {}
        autocor_matrices[sig] = {}
        gene_median_autocor[sig] = {}

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

            # Compute median autocorrelation per gene (for RankProd)
            for i, gene in enumerate(genes_present):
                if gene not in gene_median_autocor[sig]:
                    gene_median_autocor[sig][gene] = {}
                gene_median_autocor[sig][gene][ds] = float(np.nanmedian(autocors[i, :]))

            if autocors.shape[0] > 1:
                autocor_median = float(np.nanmedian(autocors))
            else:
                autocor_median = 0.0

            radar_values[sig][ds] = {"autocor_median": autocor_median}

    # RankProd analysis (only if >1 dataset)
    if compute_rank_product and len(names_datasets) > 1:
        for sig in names_sigs:
            genes = list(gene_median_autocor[sig].keys())
            if len(genes) < 2:
                continue

            # Build matrix: genes x datasets
            overall_rank_mat = np.full((len(genes), len(names_datasets)), np.nan)
            for i, gene in enumerate(genes):
                for j, ds in enumerate(names_datasets):
                    if ds in gene_median_autocor[sig][gene]:
                        overall_rank_mat[i, j] = gene_median_autocor[sig][gene][ds]

            # Run RankProd
            rp_result = _compute_rank_product(overall_rank_mat, n_permutations=n_permutations)

            # Build output table
            table = pd.DataFrame({
                "pfp_up": rp_result["pfp_up"],
                "pfp_down": rp_result["pfp_down"],
                "pval_up": rp_result["pval_up"],
                "pval_down": rp_result["pval_down"],
                "rp_up": rp_result["rp_up"],
                "rp_down": rp_result["rp_down"],
            }, index=genes)
            table = table.sort_values("rp_up")
            rank_product_tables[sig] = table

    return {
        "radar_values": radar_values,
        "autocor_matrices": autocor_matrices,
        "rank_product_tables": rank_product_tables,
        "elapsed_seconds": time.perf_counter() - _t0,
    }
