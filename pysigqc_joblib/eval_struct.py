"""Evaluate gene signature structure — subsampled biclustering for large data.

Same API and return types as pysigqc.eval_struct.compute_struct().

Key optimization: for datasets with >10k samples, subsample columns before
SpectralBiclustering to avoid O(n^3) eigendecomposition.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.cluster import SpectralBiclustering

from ._core import to_numpy, gene_indices, z_transform_matrix
from .utils import gene_intersection

BICLUST_MAX_SAMPLES = 10_000


def _compute_biclust(arr: np.ndarray, inter: list[str], col_names: list[str],
                     rng: np.random.Generator | None = None) -> dict:
    """Z-transform, binarize, and bicluster for one sig-dataset pair.

    arr: (n_genes, n_samples) submatrix for signature genes.
    """
    sig_scores = arr.copy()
    sig_scores[~np.isfinite(sig_scores)] = np.nan

    # Vectorized z-transform
    z_scores = z_transform_matrix(sig_scores)

    # Binarize at midpoint of z-score range
    z_flat = z_scores[np.isfinite(z_scores)]
    if len(z_flat) == 0:
        threshold = 0.0
    else:
        threshold = z_flat.min() + (z_flat.max() - z_flat.min()) / 2

    binarized = np.where(z_scores > threshold, 1.0, 0.0)
    binarized[np.isnan(sig_scores)] = 0.0

    # Subsample columns for large datasets
    n_samples = binarized.shape[1]
    used_cols = col_names
    if n_samples > BICLUST_MAX_SAMPLES and rng is not None:
        sample_idx = rng.choice(n_samples, size=BICLUST_MAX_SAMPLES, replace=False)
        sample_idx.sort()
        binarized_for_biclust = binarized[:, sample_idx]
        used_cols = [col_names[i] for i in sample_idx]
    else:
        binarized_for_biclust = binarized

    n_clusters = min(2, min(binarized_for_biclust.shape))
    biclust_result = None
    n_biclusters = 0
    try:
        if binarized_for_biclust.shape[0] >= 2 and binarized_for_biclust.shape[1] >= 2:
            model = SpectralBiclustering(n_clusters=n_clusters, random_state=42)
            model.fit(binarized_for_biclust)
            biclust_result = model
            n_biclusters = n_clusters
    except Exception:
        pass

    return {
        "z_scores": pd.DataFrame(z_scores, index=inter, columns=col_names),
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
    covariates: dict | None = None,
) -> dict:
    """Compute structure metrics for each signature-dataset pair.

    Args:
        covariates: Optional dict of dataset -> {'annotations': DataFrame/Series, 'colors': dict}
                   for heatmap annotations
    """
    rng = np.random.default_rng(42)

    # Pre-convert datasets
    ds_cache: dict = {}
    for ds in names_datasets:
        ds_cache[ds] = to_numpy(mRNA_expr_matrix[ds])

    # --- Step 1: Build union of gene names per signature ---
    all_row_names: dict = {}
    for sig in names_sigs:
        gene_sig = gene_sigs_list[sig]
        all_genes: list[str] = []
        for ds in names_datasets:
            _, row_names, _ = ds_cache[ds]
            inter = gene_intersection(gene_sig, mRNA_expr_matrix[ds])
            all_genes.extend(inter)
        all_row_names[sig] = list(dict.fromkeys(all_genes))

    # --- Step 2: Build padded expression matrices ---
    sig_scores_all_mats: dict = {}
    for sig in names_sigs:
        sig_scores_all_mats[sig] = {}
        for ds in names_datasets:
            arr, row_names, col_names = ds_cache[ds]
            inter = gene_intersection(gene_sigs_list[sig], mRNA_expr_matrix[ds])
            sig_idx = gene_indices(inter, row_names)
            sig_scores = pd.DataFrame(arr[sig_idx], index=inter, columns=col_names)

            # Pad missing genes with NA rows (BUG-8 fix)
            rows_needed = [g for g in all_row_names[sig] if g not in inter]
            if rows_needed:
                na_rows = pd.DataFrame(np.nan, index=rows_needed, columns=col_names)
                sig_scores = pd.concat([sig_scores, na_rows])

            sig_scores_all_mats[sig][ds] = sig_scores

    # --- Step 3: Run biclustering ---
    biclust_results: dict = {}
    any_biclusters = False
    for sig in names_sigs:
        biclust_results[sig] = {}
        for ds in names_datasets:
            arr, row_names, col_names = ds_cache[ds]
            inter = gene_intersection(gene_sigs_list[sig], mRNA_expr_matrix[ds])
            sig_idx = gene_indices(inter, row_names)
            sig_arr = arr[sig_idx]

            bc = _compute_biclust(sig_arr, inter, col_names, rng=rng)
            biclust_results[sig][ds] = bc
            if bc["n_biclusters"] > 1:
                any_biclusters = True

    return {
        "sig_scores_all_mats": sig_scores_all_mats,
        "all_row_names": all_row_names,
        "biclust_results": biclust_results,
        "any_biclusters": any_biclusters,
        "covariates": covariates,  # Pass through for plotting
    }
