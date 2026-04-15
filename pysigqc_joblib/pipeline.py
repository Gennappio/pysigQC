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
from .eval_struct import compute_struct
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
    do_negative_control: bool = False,
    num_resampling: int = 50,
    thresholds: dict[str, float] | list[float] | None = None,
    covariates: dict | None = None,
) -> dict:
    """Run the full QC pipeline.

    Parallelizes across 6 independent compute modules, then assembles
    the radar chart from their outputs. Optionally runs negative/permutation
    controls.

    Args:
        gene_sigs_list: dict of signature name -> gene list
        names_sigs: ordered list of signature names
        mRNA_expr_matrix: dict of dataset name -> DataFrame (genes x samples)
        names_datasets: ordered list of dataset names
        n_jobs: number of parallel jobs (-1 = all cores)
        verbose: if True, print progress
        do_negative_control: if True, run negative and permutation controls
        num_resampling: number of resamplings for negative control (default 50)
        thresholds: expression thresholds per dataset (dict or list, default: median)
        covariates: dict of dataset -> {annotations, colors} for struct heatmaps

    Returns dict with:
        var_result, expr_result, compact_result, stan_result, metrics_result,
        struct_result: Individual module results
        radar_result: Assembled radar chart with output_table
        radar_values: Merged per-sig per-dataset metric dicts
        negative_result: Negative/permutation control results (if requested)
    """
    _t0 = time.perf_counter()

    # Convert thresholds list to dict if needed
    if isinstance(thresholds, list):
        if len(thresholds) != len(names_datasets):
            raise ValueError(
                f"Number of thresholds ({len(thresholds)}) must match "
                f"number of datasets ({len(names_datasets)})"
            )
        thresholds = dict(zip(names_datasets, thresholds))

    # Run modules that don't need special parameters
    modules_standard = [
        ("var", compute_var),
        ("compactness", compute_compactness),
        ("stan", compute_stan),
        ("metrics", compute_metrics),
    ]

    # Run standard modules in parallel
    results_standard = Parallel(n_jobs=min(n_jobs if n_jobs > 0 else 4, 4), prefer="threads")(
        delayed(_compute_module)(
            fn, gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets
        )
        for name, fn in modules_standard
    )
    var_r, compact_r, stan_r, metrics_r = results_standard

    # Run expr with thresholds parameter
    expr_r = compute_expr(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets, thresholds=thresholds)

    # Run struct with covariates parameter
    struct_r = compute_struct(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets, covariates=covariates)

    if verbose:
        print("All compute modules complete. Assembling radar chart...")

    # Merge radar values (struct has no radar_values)
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

    result = {
        "var_result": var_r,
        "expr_result": expr_r,
        "compact_result": compact_r,
        "stan_result": stan_r,
        "metrics_result": metrics_r,
        "struct_result": struct_r,
        "radar_result": radar_result,
        "radar_values": radar_values,
        "elapsed_seconds": time.perf_counter() - _t0,
    }

    # Optionally run negative/permutation controls
    if do_negative_control:
        if verbose:
            print(f"Running negative/permutation controls ({num_resampling} resamplings)...")
        from pysigqc.negative_control import run_negative_control
        neg_result = run_negative_control(
            gene_sigs_list, mRNA_expr_matrix,
            num_resampling=num_resampling,
        )
        result["negative_result"] = neg_result
        result["elapsed_seconds"] = time.perf_counter() - _t0

    return result
