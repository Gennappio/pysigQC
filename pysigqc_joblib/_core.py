"""Vectorized numpy primitives for hot paths — no pandas in compute loops.

All functions operate on raw numpy arrays. DataFrame conversion happens once
at module entry via to_numpy().
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as sp_stats


def to_numpy(df: pd.DataFrame) -> tuple[np.ndarray, list[str], list[str]]:
    """Convert DataFrame to (values, row_names, col_names) once."""
    return df.values.astype(np.float64, copy=False), list(df.index), list(df.columns)


def gene_indices(signature: list[str], row_names: list[str]) -> np.ndarray:
    """Return integer indices of signature genes present in row_names."""
    name_to_idx = {name: i for i, name in enumerate(row_names)}
    return np.array([name_to_idx[g] for g in signature if g in name_to_idx], dtype=np.intp)


def nanstd_rows(arr: np.ndarray, ddof: int = 1) -> np.ndarray:
    """Row-wise std with ddof, pure numpy. All-NaN rows get NaN."""
    with np.errstate(invalid="ignore"):
        return np.nanstd(arr, axis=1, ddof=ddof)


def nanmean_rows(arr: np.ndarray) -> np.ndarray:
    """Row-wise mean, pure numpy."""
    with np.errstate(invalid="ignore"):
        return np.nanmean(arr, axis=1)


def nanmedian_cols(arr: np.ndarray) -> np.ndarray:
    """Column-wise median (across genes/rows for each sample/column)."""
    with np.errstate(invalid="ignore"):
        return np.nanmedian(arr, axis=0)


def nanmedian_rows(arr: np.ndarray) -> np.ndarray:
    """Row-wise median."""
    with np.errstate(invalid="ignore"):
        return np.nanmedian(arr, axis=1)


def z_transform_matrix(arr: np.ndarray) -> np.ndarray:
    """Vectorized z-transform all rows at once, zero-variance rows become all-zero."""
    with np.errstate(invalid="ignore"):
        means = np.nanmean(arr, axis=1, keepdims=True)
        stds = np.nanstd(arr, axis=1, ddof=1, keepdims=True)
    # Guard: zero-variance or all-NaN rows
    bad = (stds == 0) | np.isnan(stds)
    stds_safe = np.where(bad, 1.0, stds)
    result = (arr - means) / stds_safe
    # Zero out rows that had zero variance (instead of 0/1 = 0, handles NaN case)
    result[bad.squeeze(axis=1)] = 0.0
    return result


def spearman_matrix(arr: np.ndarray) -> np.ndarray:
    """Spearman correlation matrix for rows of arr.

    Fast path: if no NaN, rank all rows at once and use np.corrcoef.
    Slow path: pairwise spearmanr for rows containing NaN.

    arr shape: (n_genes, n_samples)
    Returns: (n_genes, n_genes) correlation matrix
    """
    n = arr.shape[0]
    if n <= 1:
        return np.ones((n, n))

    has_nan = np.isnan(arr).any()

    if not has_nan:
        # Fast path: rank once, corrcoef once
        ranks = sp_stats.rankdata(arr, axis=1).astype(np.float64)
        with np.errstate(invalid="ignore"):
            corr = np.corrcoef(ranks)
        # Constant rows produce NaN on diagonal — force diagonal to 1.0
        # to match the pairwise spearmanr behavior (identity on diagonal)
        np.fill_diagonal(corr, 1.0)
        return corr

    # Slow path: pairwise, handling NaN per pair
    corr = np.eye(n, dtype=np.float64)
    for i in range(n):
        for j in range(i + 1, n):
            # Use only positions where both genes are non-NaN
            mask = ~(np.isnan(arr[i]) | np.isnan(arr[j]))
            if mask.sum() < 3:
                corr[i, j] = corr[j, i] = np.nan
            else:
                rho, _ = sp_stats.spearmanr(arr[i, mask], arr[j, mask])
                corr[i, j] = corr[j, i] = rho
    return corr


def skew_rows(arr: np.ndarray, bias: bool = True) -> np.ndarray:
    """Row-wise skewness, matching scipy.stats.skew with nan_policy='omit'."""
    return sp_stats.skew(arr, axis=1, bias=bias, nan_policy="omit")


def rows_without_nan(arr: np.ndarray) -> np.ndarray:
    """Boolean mask of rows that have no NaN values."""
    return ~np.isnan(arr).any(axis=1)
