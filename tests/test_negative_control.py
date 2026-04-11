"""Tests for pysigqc.negative_control — port of R sigsQcNegativeControl."""

import numpy as np
import pandas as pd
import pytest

from pysigqc.negative_control import _compute_qc_metrics, run_negative_control


def test_compute_qc_metrics_structure(signatures, datasets, names_sigs, names_datasets):
    result = _compute_qc_metrics(signatures, names_sigs, datasets, names_datasets)
    assert "radar_values" in result
    assert "output_table" in result
    assert result["output_table"].shape == (4, 14)  # 2 sigs * 2 datasets


def test_compute_qc_metrics_non_negative(signatures, datasets, names_sigs, names_datasets):
    result = _compute_qc_metrics(signatures, names_sigs, datasets, names_datasets)
    assert (result["output_table"].values >= 0).all()


def test_run_negative_control_small(signatures, datasets):
    """Run with very few resamplings to check structure."""
    result = run_negative_control(
        {"compact_sig": signatures["compact_sig"]},
        {"dataset_A": datasets["dataset_A"]},
        num_resampling=3,
        seed=42,
    )

    assert "negative_controls" in result
    assert "permutation_controls" in result

    nc = result["negative_controls"]["dataset_A"]["compact_sig"]
    assert "summary" in nc
    assert "metrics_table" in nc
    assert nc["metrics_table"].shape[0] == 3  # 3 resamplings
    assert nc["metrics_table"].shape[1] == 14  # 14 metrics
    assert nc["summary"].shape == (6, 14)  # mean + 5 quantiles

    pc = result["permutation_controls"]["dataset_A"]["compact_sig"]
    assert pc["metrics_table"].shape[0] == 3
    assert pc["summary"].shape == (6, 14)


def test_negative_control_deterministic(signatures, datasets):
    """Same seed should produce same results."""
    r1 = run_negative_control(
        {"compact_sig": signatures["compact_sig"]},
        {"dataset_A": datasets["dataset_A"]},
        num_resampling=3,
        seed=42,
    )
    r2 = run_negative_control(
        {"compact_sig": signatures["compact_sig"]},
        {"dataset_A": datasets["dataset_A"]},
        num_resampling=3,
        seed=42,
    )

    pd.testing.assert_frame_equal(
        r1["negative_controls"]["dataset_A"]["compact_sig"]["metrics_table"],
        r2["negative_controls"]["dataset_A"]["compact_sig"]["metrics_table"],
    )
