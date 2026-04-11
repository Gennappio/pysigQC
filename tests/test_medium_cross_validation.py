"""Cross-validation tests on the medium-scale fixture (100 genes, 50 samples, 3 datasets, 5 sigs).

Tests numerical stability at more realistic dimensions and exercises:
- batch effects between datasets
- zero-variance and near-zero-variance genes
- scattered, block, and full-row NA patterns
- genes missing from one dataset but not others
- signature gene overlap
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pysigqc.eval_var import compute_var
from pysigqc.eval_expr import compute_expr
from pysigqc.eval_compactness import compute_compactness
from pysigqc.eval_stan import compute_stan
from pysigqc.compare_metrics import compute_metrics
from pysigqc.radar_chart import compute_radar

FIXTURES_DIR = Path(__file__).resolve().parent.parent / "sigQC-master" / "tests" / "fixtures"
REF_DIR = FIXTURES_DIR / "reference_outputs_medium"


@pytest.fixture(scope="module")
def medium_datasets() -> dict[str, pd.DataFrame]:
    ds = {}
    for name in ["dataset_A", "dataset_B", "dataset_C"]:
        f = FIXTURES_DIR / f"fixture_medium_{name}.csv"
        if not f.exists():
            pytest.skip(f"Medium fixture not found: {f}")
        ds[name] = pd.read_csv(f, index_col=0)
    return ds


@pytest.fixture(scope="module")
def medium_signatures() -> dict[str, list[str]]:
    f = FIXTURES_DIR / "fixture_medium_signatures.csv"
    if not f.exists():
        pytest.skip("Medium signatures file not found")
    df = pd.read_csv(f)
    sigs: dict[str, list[str]] = {}
    for sig_name, group in df.groupby("signature"):
        sigs[sig_name] = group["gene"].tolist()
    return sigs


@pytest.fixture(scope="module")
def medium_names_sigs() -> list[str]:
    return ["tight_coexpr", "loose_coexpr", "mixed_quality", "high_var", "random_baseline"]


@pytest.fixture(scope="module")
def medium_names_datasets() -> list[str]:
    return ["dataset_A", "dataset_B", "dataset_C"]


# ---------------------------------------------------------------------------
# Per-module cross-validation
# ---------------------------------------------------------------------------

def _check_radar_against_ref(py_vals: dict, ref_dir: Path, prefix: str,
                              sig: str, ds: str, rtol: float = 1e-3):
    ref_file = ref_dir / f"{prefix}_radar_{sig}_{ds}.csv"
    if not ref_file.exists():
        pytest.skip(f"No ref: {ref_file}")
    ref = pd.read_csv(ref_file)
    for col in ref.columns:
        r_val = ref[col].iloc[0]
        py_val = py_vals.get(col, 0.0)
        if np.isnan(r_val) and np.isnan(py_val):
            continue
        # PCA sign ambiguity: compare abs for rho involving PCA
        if "pca1" in col:
            np.testing.assert_allclose(
                abs(py_val), abs(r_val), rtol=rtol, atol=1e-8,
                err_msg=f"|{col}| mismatch for {sig}/{ds}: Py={py_val}, R={r_val}"
            )
        else:
            np.testing.assert_allclose(
                py_val, r_val, rtol=rtol, atol=1e-8,
                err_msg=f"{col} mismatch for {sig}/{ds}: Py={py_val}, R={r_val}"
            )


def test_var_medium(medium_signatures, medium_datasets, medium_names_sigs, medium_names_datasets):
    result = compute_var(medium_signatures, medium_names_sigs, medium_datasets, medium_names_datasets)
    for sig in medium_names_sigs:
        for ds in medium_names_datasets:
            rv = result["radar_values"][sig][ds]
            assert all(0 <= v <= 1 for v in rv.values() if not np.isnan(v))
            _check_radar_against_ref(rv, REF_DIR, "var", sig, ds)


def test_expr_medium(medium_signatures, medium_datasets, medium_names_sigs, medium_names_datasets):
    result = compute_expr(medium_signatures, medium_names_sigs, medium_datasets, medium_names_datasets)
    for sig in medium_names_sigs:
        for ds in medium_names_datasets:
            rv = result["radar_values"][sig][ds]
            assert all(0 <= v <= 1 for v in rv.values() if not np.isnan(v))
            _check_radar_against_ref(rv, REF_DIR, "expr", sig, ds)


def test_compactness_medium(medium_signatures, medium_datasets, medium_names_sigs, medium_names_datasets):
    result = compute_compactness(medium_signatures, medium_names_sigs, medium_datasets, medium_names_datasets)
    for sig in medium_names_sigs:
        for ds in medium_names_datasets:
            rv = result["radar_values"][sig][ds]
            _check_radar_against_ref(rv, REF_DIR, "compact", sig, ds)


def test_stan_medium(medium_signatures, medium_datasets, medium_names_sigs, medium_names_datasets):
    result = compute_stan(medium_signatures, medium_names_sigs, medium_datasets, medium_names_datasets)
    for sig in medium_names_sigs:
        for ds in medium_names_datasets:
            rv = result["radar_values"][sig][ds]
            _check_radar_against_ref(rv, REF_DIR, "stan", sig, ds)


def test_metrics_medium(medium_signatures, medium_datasets, medium_names_sigs, medium_names_datasets):
    result = compute_metrics(medium_signatures, medium_names_sigs, medium_datasets, medium_names_datasets)
    for sig in medium_names_sigs:
        for ds in medium_names_datasets:
            rv = result["radar_values"][sig][ds]
            _check_radar_against_ref(rv, REF_DIR, "metrics", sig, ds, rtol=0.05)


def test_full_radar_medium(medium_signatures, medium_datasets, medium_names_sigs, medium_names_datasets):
    """Full pipeline: all 14 metrics assembled into radar chart, compared to R."""
    var_r = compute_var(medium_signatures, medium_names_sigs, medium_datasets, medium_names_datasets)
    expr_r = compute_expr(medium_signatures, medium_names_sigs, medium_datasets, medium_names_datasets)
    compact_r = compute_compactness(medium_signatures, medium_names_sigs, medium_datasets, medium_names_datasets)
    stan_r = compute_stan(medium_signatures, medium_names_sigs, medium_datasets, medium_names_datasets)
    metrics_r = compute_metrics(medium_signatures, medium_names_sigs, medium_datasets, medium_names_datasets)

    radar_values = {}
    for sig in medium_names_sigs:
        radar_values[sig] = {}
        for ds in medium_names_datasets:
            vals = {}
            vals.update(var_r["radar_values"][sig][ds])
            vals.update(expr_r["radar_values"][sig][ds])
            vals.update(compact_r["radar_values"][sig][ds])
            vals.update(metrics_r["radar_values"][sig][ds])
            vals.update(stan_r["radar_values"][sig][ds])
            radar_values[sig][ds] = vals

    result = compute_radar(radar_values, medium_names_sigs, medium_names_datasets)

    ref_file = REF_DIR / "radar_output_table.csv"
    if not ref_file.exists():
        pytest.skip("No medium radar reference")
    ref = pd.read_csv(ref_file, index_col=0)

    assert result["output_table"].shape == ref.shape

    for col in ref.columns:
        for idx in ref.index:
            if idx not in result["output_table"].index:
                continue
            r_val = ref.loc[idx, col]
            py_val = result["output_table"].loc[idx, col]
            if np.isnan(r_val) and np.isnan(py_val):
                continue
            if "pca1" in col:
                np.testing.assert_allclose(abs(py_val), abs(r_val), rtol=0.05, atol=1e-8,
                                           err_msg=f"|{col}| at {idx}")
            else:
                np.testing.assert_allclose(py_val, r_val, rtol=1e-3, atol=1e-8,
                                           err_msg=f"{col} at {idx}")
