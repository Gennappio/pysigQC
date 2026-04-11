"""Full QC pipeline orchestrator with parallel execution.

Parallelizes across (signature, dataset) pairs using joblib.
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from .eval_var import compute_var
from .eval_expr import compute_expr
from .eval_compactness import compute_compactness
from .eval_stan import compute_stan
from .compare_metrics import compute_metrics
from .radar_chart import compute_radar


def _compute_module(module_fn, gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets):
    """Run a single compute module. Used as a joblib target."""
    return module_fn(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)


def run_pipeline(
    gene_sigs_list: dict[str, list[str]],
    names_sigs: list[str],
    mRNA_expr_matrix: dict[str, pd.DataFrame],
    names_datasets: list[str],
    n_jobs: int = -1,
    verbose: bool = False,
) -> dict:
    """Run the full QC pipeline.

    Parallelizes across the 5 independent compute modules, then assembles
    the radar chart from their outputs.

    Args:
        gene_sigs_list: dict of signature name -> gene list
        names_sigs: ordered list of signature names
        mRNA_expr_matrix: dict of dataset name -> DataFrame (genes x samples)
        names_datasets: ordered list of dataset names
        n_jobs: number of parallel jobs (-1 = all cores)
        verbose: if True, print progress

    Returns dict with:
        var_result, expr_result, compact_result, stan_result, metrics_result:
            Individual module results
        radar_result: Assembled radar chart with output_table
        radar_values: Merged per-sig per-dataset metric dicts
    """
    _t0 = time.perf_counter()

    modules = [
        ("var", compute_var),
        ("expr", compute_expr),
        ("compactness", compute_compactness),
        ("stan", compute_stan),
        ("metrics", compute_metrics),
    ]

    # Run all 5 modules in parallel
    results = Parallel(n_jobs=min(n_jobs if n_jobs > 0 else 5, 5), prefer="threads")(
        delayed(_compute_module)(
            fn, gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets
        )
        for name, fn in modules
    )

    var_r, expr_r, compact_r, stan_r, metrics_r = results

    if verbose:
        print("All compute modules complete. Assembling radar chart...")

    # Merge radar values
    radar_values: dict = {}
    for sig in names_sigs:
        radar_values[sig] = {}
        for ds in names_datasets:
            vals: dict = {}
            vals.update(var_r["radar_values"][sig][ds])
            vals.update(expr_r["radar_values"][sig][ds])
            vals.update(compact_r["radar_values"][sig][ds])
            vals.update(metrics_r["radar_values"][sig][ds])
            vals.update(stan_r["radar_values"][sig][ds])
            radar_values[sig][ds] = vals

    radar_result = compute_radar(radar_values, names_sigs, names_datasets)

    return {
        "var_result": var_r,
        "expr_result": expr_r,
        "compact_result": compact_r,
        "stan_result": stan_r,
        "metrics_result": metrics_r,
        "radar_result": radar_result,
        "radar_values": radar_values,
        "elapsed_seconds": time.perf_counter() - _t0,
    }
