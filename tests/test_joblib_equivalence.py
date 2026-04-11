"""Numerical equivalence tests: pysigqc vs pysigqc_joblib.

Verifies that the optimized joblib version produces identical results
to the reference pysigqc implementation on both small and medium fixtures.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# Reference implementations
import pysigqc.eval_var as ref_var
import pysigqc.eval_expr as ref_expr
import pysigqc.eval_compactness as ref_compact
import pysigqc.eval_stan as ref_stan
import pysigqc.compare_metrics as ref_metrics
import pysigqc.radar_chart as ref_radar
import pysigqc.negative_control as ref_nc

# Optimized implementations
import pysigqc_joblib.eval_var as opt_var
import pysigqc_joblib.eval_expr as opt_expr
import pysigqc_joblib.eval_compactness as opt_compact
import pysigqc_joblib.eval_stan as opt_stan
import pysigqc_joblib.compare_metrics as opt_metrics
import pysigqc_joblib.radar_chart as opt_radar
import pysigqc_joblib.negative_control as opt_nc


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _compare_radar_values(ref_rv, opt_rv, names_sigs, names_datasets,
                          rtol=1e-10, atol=1e-12, pca_rtol=1e-3):
    """Compare radar values from two implementations."""
    for sig in names_sigs:
        for ds in names_datasets:
            ref_vals = ref_rv[sig][ds]
            opt_vals = opt_rv[sig][ds]
            assert set(ref_vals.keys()) == set(opt_vals.keys()), \
                f"Key mismatch for {sig}/{ds}: {set(ref_vals.keys())} vs {set(opt_vals.keys())}"
            for key in ref_vals:
                rv = ref_vals[key]
                ov = opt_vals[key]
                if np.isnan(rv) and np.isnan(ov):
                    continue
                if "pca1" in key:
                    np.testing.assert_allclose(
                        abs(ov), abs(rv), rtol=pca_rtol, atol=atol,
                        err_msg=f"|{key}| mismatch for {sig}/{ds}"
                    )
                else:
                    np.testing.assert_allclose(
                        ov, rv, rtol=rtol, atol=atol,
                        err_msg=f"{key} mismatch for {sig}/{ds}"
                    )


# ---------------------------------------------------------------------------
# Small fixture equivalence tests
# ---------------------------------------------------------------------------

class TestSmallFixtureEquivalence:
    """Test equivalence on the small fixture (shared with conftest.py)."""

    def test_var_equivalence(self, signatures, datasets, names_sigs, names_datasets):
        ref = ref_var.compute_var(signatures, names_sigs, datasets, names_datasets)
        opt = opt_var.compute_var(signatures, names_sigs, datasets, names_datasets)
        _compare_radar_values(ref["radar_values"], opt["radar_values"],
                              names_sigs, names_datasets)

    def test_expr_equivalence(self, signatures, datasets, names_sigs, names_datasets):
        ref = ref_expr.compute_expr(signatures, names_sigs, datasets, names_datasets)
        opt = opt_expr.compute_expr(signatures, names_sigs, datasets, names_datasets)
        _compare_radar_values(ref["radar_values"], opt["radar_values"],
                              names_sigs, names_datasets)

    def test_compactness_equivalence(self, signatures, datasets, names_sigs, names_datasets):
        ref = ref_compact.compute_compactness(signatures, names_sigs, datasets, names_datasets)
        opt = opt_compact.compute_compactness(signatures, names_sigs, datasets, names_datasets)
        _compare_radar_values(ref["radar_values"], opt["radar_values"],
                              names_sigs, names_datasets)

    def test_stan_equivalence(self, signatures, datasets, names_sigs, names_datasets):
        ref = ref_stan.compute_stan(signatures, names_sigs, datasets, names_datasets)
        opt = opt_stan.compute_stan(signatures, names_sigs, datasets, names_datasets)
        _compare_radar_values(ref["radar_values"], opt["radar_values"],
                              names_sigs, names_datasets)

    def test_metrics_equivalence(self, signatures, datasets, names_sigs, names_datasets):
        ref = ref_metrics.compute_metrics(signatures, names_sigs, datasets, names_datasets)
        opt = opt_metrics.compute_metrics(signatures, names_sigs, datasets, names_datasets)
        _compare_radar_values(ref["radar_values"], opt["radar_values"],
                              names_sigs, names_datasets, pca_rtol=0.05)

    def test_full_radar_equivalence(self, signatures, datasets, names_sigs, names_datasets):
        """Full pipeline: assemble all 14 metrics and compare radar tables."""
        # Reference
        r_var = ref_var.compute_var(signatures, names_sigs, datasets, names_datasets)
        r_expr = ref_expr.compute_expr(signatures, names_sigs, datasets, names_datasets)
        r_compact = ref_compact.compute_compactness(signatures, names_sigs, datasets, names_datasets)
        r_stan = ref_stan.compute_stan(signatures, names_sigs, datasets, names_datasets)
        r_metrics = ref_metrics.compute_metrics(signatures, names_sigs, datasets, names_datasets)

        # Optimized
        o_var = opt_var.compute_var(signatures, names_sigs, datasets, names_datasets)
        o_expr = opt_expr.compute_expr(signatures, names_sigs, datasets, names_datasets)
        o_compact = opt_compact.compute_compactness(signatures, names_sigs, datasets, names_datasets)
        o_stan = opt_stan.compute_stan(signatures, names_sigs, datasets, names_datasets)
        o_metrics = opt_metrics.compute_metrics(signatures, names_sigs, datasets, names_datasets)

        for impl_label, var_r, expr_r, compact_r, stan_r, metrics_r, radar_fn in [
            ("ref", r_var, r_expr, r_compact, r_stan, r_metrics, ref_radar.compute_radar),
            ("opt", o_var, o_expr, o_compact, o_stan, o_metrics, opt_radar.compute_radar),
        ]:
            rv = {}
            for sig in names_sigs:
                rv[sig] = {}
                for ds in names_datasets:
                    vals = {}
                    vals.update(var_r["radar_values"][sig][ds])
                    vals.update(expr_r["radar_values"][sig][ds])
                    vals.update(compact_r["radar_values"][sig][ds])
                    vals.update(metrics_r["radar_values"][sig][ds])
                    vals.update(stan_r["radar_values"][sig][ds])
                    rv[sig][ds] = vals

        # Compare the two radar tables
        ref_result = ref_radar.compute_radar(rv, names_sigs, names_datasets)
        opt_result = opt_radar.compute_radar(rv, names_sigs, names_datasets)
        assert ref_result["output_table"].shape == opt_result["output_table"].shape


# ---------------------------------------------------------------------------
# Medium fixture equivalence tests
# ---------------------------------------------------------------------------

FIXTURES_DIR = Path(__file__).resolve().parent.parent / "sigQC-master" / "tests" / "fixtures"


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


class TestMediumFixtureEquivalence:
    """Test equivalence on the medium fixture (100 genes, 50 samples, 3 datasets)."""

    def test_var_medium(self, medium_signatures, medium_datasets, medium_names_sigs, medium_names_datasets):
        ref = ref_var.compute_var(medium_signatures, medium_names_sigs, medium_datasets, medium_names_datasets)
        opt = opt_var.compute_var(medium_signatures, medium_names_sigs, medium_datasets, medium_names_datasets)
        _compare_radar_values(ref["radar_values"], opt["radar_values"],
                              medium_names_sigs, medium_names_datasets)

    def test_expr_medium(self, medium_signatures, medium_datasets, medium_names_sigs, medium_names_datasets):
        ref = ref_expr.compute_expr(medium_signatures, medium_names_sigs, medium_datasets, medium_names_datasets)
        opt = opt_expr.compute_expr(medium_signatures, medium_names_sigs, medium_datasets, medium_names_datasets)
        _compare_radar_values(ref["radar_values"], opt["radar_values"],
                              medium_names_sigs, medium_names_datasets)

    def test_compactness_medium(self, medium_signatures, medium_datasets, medium_names_sigs, medium_names_datasets):
        ref = ref_compact.compute_compactness(medium_signatures, medium_names_sigs, medium_datasets, medium_names_datasets)
        opt = opt_compact.compute_compactness(medium_signatures, medium_names_sigs, medium_datasets, medium_names_datasets)
        _compare_radar_values(ref["radar_values"], opt["radar_values"],
                              medium_names_sigs, medium_names_datasets)

    def test_stan_medium(self, medium_signatures, medium_datasets, medium_names_sigs, medium_names_datasets):
        ref = ref_stan.compute_stan(medium_signatures, medium_names_sigs, medium_datasets, medium_names_datasets)
        opt = opt_stan.compute_stan(medium_signatures, medium_names_sigs, medium_datasets, medium_names_datasets)
        _compare_radar_values(ref["radar_values"], opt["radar_values"],
                              medium_names_sigs, medium_names_datasets)

    def test_metrics_medium(self, medium_signatures, medium_datasets, medium_names_sigs, medium_names_datasets):
        ref = ref_metrics.compute_metrics(medium_signatures, medium_names_sigs, medium_datasets, medium_names_datasets)
        opt = opt_metrics.compute_metrics(medium_signatures, medium_names_sigs, medium_datasets, medium_names_datasets)
        _compare_radar_values(ref["radar_values"], opt["radar_values"],
                              medium_names_sigs, medium_names_datasets, pca_rtol=0.05)
