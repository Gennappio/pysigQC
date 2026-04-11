"""Evaluate compactness (internal coherence) of gene signatures via autocorrelation.

Port of R_refactored/eval_compactness_loc.R — compute_compactness() only.
Produces 1 radar metric: autocor_median (median of Spearman gene-gene correlation matrix).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as sp_stats

from .utils import gene_intersection


def compute_compactness(
    gene_sigs_list: dict[str, list[str]],
    names_sigs: list[str],
    mRNA_expr_matrix: dict[str, pd.DataFrame],
    names_datasets: list[str],
) -> dict:
    """Compute compactness metrics for each signature-dataset pair.

    Returns dict with keys:
        radar_values: nested dict [sig][dataset] -> {"autocor_median": val}
        autocor_matrices: nested dict [sig][dataset] -> gene-gene Spearman correlation matrix
        rank_product_tables: dict (empty in Python — RankProd not ported)
    """
    radar_values: dict = {}
    autocor_matrices: dict = {}

    for sig in names_sigs:
        gene_sig = gene_sigs_list[sig]
        radar_values[sig] = {}
        autocor_matrices[sig] = {}

        for ds in names_datasets:
            data_matrix = mRNA_expr_matrix[ds]
            inter = gene_intersection(gene_sig, data_matrix)

            # Get sig gene data, drop genes with any NA (matches R's na.omit)
            sig_df = data_matrix.loc[inter].dropna(axis=0, how="any")
            genes_present = list(sig_df.index)
            sig_data = sig_df.values.astype(float)

            # Spearman correlation between genes (rows)
            # R does: cor(t(na.omit(data.matrix[inter,])), method='spearman')
            # which computes pairwise correlations — zero-variance genes get NaN
            # for their pairs but don't affect other pairs.
            n_genes = sig_data.shape[0]
            if n_genes > 1:
                # Build pairwise correlation matrix to match R's cor() behavior
                autocors = np.eye(n_genes)
                for gi in range(n_genes):
                    for gj in range(gi + 1, n_genes):
                        rho, _ = sp_stats.spearmanr(sig_data[gi], sig_data[gj])
                        autocors[gi, gj] = rho
                        autocors[gj, gi] = rho
            else:
                autocors = np.array([[1.0]])

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
    }
