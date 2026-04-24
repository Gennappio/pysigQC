"""Evaluate gene signature structure via hierarchical clustering and biclustering.

Port of R_refactored/eval_struct_loc.R — compute_struct() only.
Uses scipy for hierarchical clustering and sklearn for spectral biclustering
as a substitute for R's biclust BCCC method.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from sklearn.cluster import SpectralBiclustering

from .utils import gene_intersection, z_transform


def _compute_biclust(data_matrix: pd.DataFrame, inter: list[str]) -> dict:
    """Z-transform, binarize, and bicluster for one sig-dataset pair.

    Mirrors R's .compute_biclust() helper.
    """
    sig_scores = data_matrix.loc[inter].copy().values.astype(float)
    sig_scores[~np.isfinite(sig_scores)] = np.nan

    # Z-transform each gene (with zero-variance guard)
    for i, gene in enumerate(inter):
        sig_scores[i] = z_transform(sig_scores[i])

    # Binarize at midpoint of z-score range
    z_flat = sig_scores[np.isfinite(sig_scores)]
    if len(z_flat) == 0:
        threshold = 0.0
    else:
        threshold = z_flat.min() + (z_flat.max() - z_flat.min()) / 2

    binarized = np.where(sig_scores > threshold, 1, 0).astype(float)
    binarized[np.isnan(sig_scores)] = 0  # fill NaN with 0 for biclustering

    # Biclustering (using SpectralBiclustering as Python analog of BCCC)
    n_clusters = min(2, min(binarized.shape))
    biclust_result = None
    n_biclusters = 0
    try:
        if binarized.shape[0] >= 2 and binarized.shape[1] >= 2:
            model = SpectralBiclustering(n_clusters=n_clusters, random_state=42)
            model.fit(binarized)
            biclust_result = model
            n_biclusters = n_clusters
    except (ValueError, np.linalg.LinAlgError) as e:
        # SpectralBiclustering can fail on degenerate (low-rank / constant) data.
        # Matches R's biclust::BCCC tryCatch in R_refactored/eval_struct_loc.R.
        warnings.warn(
            f"SpectralBiclustering failed: {type(e).__name__}: {e}",
            RuntimeWarning, stacklevel=2,
        )

    return {
        "z_scores": pd.DataFrame(sig_scores, index=inter, columns=data_matrix.columns),
        "binarized": binarized,
        "biclust_result": biclust_result,
        "threshold": threshold,
        "n_biclusters": n_biclusters,
    }


def compute_struct(
    gene_sigs_list: dict[str, list[str]],
    names_sigs: list[str],
    mRNA_expr_matrix: dict[str, pd.DataFrame],
    names_datasets: list[str],
) -> dict:
    """Compute structure metrics for each signature-dataset pair.

    Returns dict with keys:
        sig_scores_all_mats: nested dict [sig][dataset] -> DataFrame (padded expression)
        all_row_names: dict [sig] -> list of union gene names across datasets
        biclust_results: nested dict [sig][dataset] -> biclust result dict
        any_biclusters: bool
    """
    # --- Step 1: Build union of gene names per signature ---
    all_row_names: dict = {}
    for sig in names_sigs:
        gene_sig = gene_sigs_list[sig]
        all_genes: list[str] = []
        for ds in names_datasets:
            data_matrix = mRNA_expr_matrix[ds]
            inter = gene_intersection(gene_sig, data_matrix)
            all_genes.extend(inter)
        all_row_names[sig] = list(dict.fromkeys(all_genes))  # unique, order preserved

    # --- Step 2: Build padded expression matrices ---
    sig_scores_all_mats: dict = {}
    for sig in names_sigs:
        sig_scores_all_mats[sig] = {}
        gene_sig = gene_sigs_list[sig]

        for ds in names_datasets:
            data_matrix = mRNA_expr_matrix[ds]
            inter = gene_intersection(gene_sig, data_matrix)
            sig_scores = data_matrix.loc[inter].copy()

            # Pad missing genes with NA rows (BUG-8 fix)
            rows_needed = [g for g in all_row_names[sig] if g not in inter]
            if rows_needed:
                na_rows = pd.DataFrame(
                    np.nan,
                    index=rows_needed,
                    columns=sig_scores.columns,
                )
                sig_scores = pd.concat([sig_scores, na_rows])

            sig_scores_all_mats[sig][ds] = sig_scores

    # --- Step 3: Run biclustering once per sig-dataset pair ---
    biclust_results: dict = {}
    any_biclusters = False
    for sig in names_sigs:
        gene_sig = gene_sigs_list[sig]
        biclust_results[sig] = {}
        for ds in names_datasets:
            data_matrix = mRNA_expr_matrix[ds]
            inter = gene_intersection(gene_sig, data_matrix)
            bc = _compute_biclust(data_matrix, inter)
            biclust_results[sig][ds] = bc
            if bc["n_biclusters"] > 1:
                any_biclusters = True

    return {
        "sig_scores_all_mats": sig_scores_all_mats,
        "all_row_names": all_row_names,
        "biclust_results": biclust_results,
        "any_biclusters": any_biclusters,
    }
