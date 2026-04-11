"""Tests for pysigqc.eval_struct — port of R compute_struct()."""

import numpy as np
import pandas as pd
import pytest

from pysigqc.eval_struct import compute_struct


def test_structure(signatures, datasets, names_sigs, names_datasets):
    result = compute_struct(signatures, names_sigs, datasets, names_datasets)
    assert set(result.keys()) == {"sig_scores_all_mats", "all_row_names",
                                   "biclust_results", "any_biclusters"}


def test_all_row_names_union(signatures, datasets, names_sigs, names_datasets):
    result = compute_struct(signatures, names_sigs, datasets, names_datasets)
    for sig in names_sigs:
        gene_sig = signatures[sig]
        union_genes = result["all_row_names"][sig]
        # Every gene in the union should be in at least one dataset
        for g in union_genes:
            found = any(g in datasets[ds].index for ds in names_datasets)
            assert found, f"Gene {g} not in any dataset for {sig}"


def test_padded_matrices_consistent_rows(signatures, datasets, names_sigs, names_datasets):
    result = compute_struct(signatures, names_sigs, datasets, names_datasets)
    for sig in names_sigs:
        expected_nrow = len(result["all_row_names"][sig])
        for ds in names_datasets:
            mat = result["sig_scores_all_mats"][sig][ds]
            assert mat.shape[0] == expected_nrow, f"Row count mismatch for {sig}/{ds}"


def test_padded_with_na_for_missing(signatures, datasets):
    """Asymmetric datasets: gene present in one but not the other gets NA padding."""
    mat_X = pd.DataFrame(
        np.random.randn(3, 10),
        index=["gene_a", "gene_b", "gene_d"],
        columns=[f"s{i}" for i in range(10)],
    )
    mat_Y = pd.DataFrame(
        np.random.randn(3, 10),
        index=["gene_b", "gene_c", "gene_d"],
        columns=[f"s{i}" for i in range(10)],
    )
    test_sigs = {"test_sig": ["gene_a", "gene_b", "gene_c"]}
    test_expr = {"dataset_X": mat_X, "dataset_Y": mat_Y}

    result = compute_struct(test_sigs, ["test_sig"], test_expr, ["dataset_X", "dataset_Y"])

    assert len(result["all_row_names"]["test_sig"]) == 3

    # dataset_X missing gene_c → NA
    mat_padded_X = result["sig_scores_all_mats"]["test_sig"]["dataset_X"]
    assert "gene_c" in mat_padded_X.index
    assert mat_padded_X.loc["gene_c"].isna().all()

    # dataset_Y missing gene_a → NA
    mat_padded_Y = result["sig_scores_all_mats"]["test_sig"]["dataset_Y"]
    assert "gene_a" in mat_padded_Y.index
    assert mat_padded_Y.loc["gene_a"].isna().all()


def test_biclust_results_populated(signatures, datasets, names_sigs, names_datasets):
    result = compute_struct(signatures, names_sigs, datasets, names_datasets)
    for sig in names_sigs:
        for ds in names_datasets:
            bc = result["biclust_results"][sig][ds]
            assert "z_scores" in bc
            assert "binarized" in bc
            assert "biclust_result" in bc
            assert "threshold" in bc


def test_any_biclusters_is_bool(signatures, datasets, names_sigs, names_datasets):
    result = compute_struct(signatures, names_sigs, datasets, names_datasets)
    assert isinstance(result["any_biclusters"], bool)
