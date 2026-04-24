"""Negative and permutation controls — fully parallelized with joblib.

Same API and return types as pysigqc.negative_control.run_negative_control().

Key optimizations:
- joblib.Parallel for independent resamplings (n_jobs cores)
- np.random.SeedSequence.spawn() for reproducible parallel RNG
- Each worker is self-contained (no shared mutable state)
- Vectorized quantile computation
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from .eval_var import compute_var
from .eval_expr import compute_expr
from .eval_compactness import compute_compactness
from .eval_stan import compute_stan
from .compare_metrics import compute_metrics
from .radar_chart import compute_radar, ALL_METRICS


def _compute_qc_metrics(
    gene_sigs_list: dict[str, list[str]],
    names_sigs: list[str],
    mRNA_expr_matrix: dict[str, pd.DataFrame],
    names_datasets: list[str],
) -> dict:
    """Compute all radar metrics in memory (no file I/O)."""
    radar_values: dict = {}
    for sig in names_sigs:
        radar_values[sig] = {}
        for ds in names_datasets:
            radar_values[sig][ds] = {}

    modules = [
        ("var", lambda: compute_var(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)),
        ("expr", lambda: compute_expr(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)),
        ("compactness", lambda: compute_compactness(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)),
        # GSVA / ssGSEA / PLAGE and Gaussian mixture models are skipped: neither
        # output is consumed by the negative/permutation control summarization,
        # which reads only the 4 radar metrics from output_table.
        ("metrics", lambda: compute_metrics(gene_sigs_list, names_sigs, mRNA_expr_matrix,
                                             names_datasets,
                                             compute_enrichment=False,
                                             fit_mixture=False)),
        ("stan", lambda: compute_stan(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)),
    ]

    # Errors from any module propagate with a stack trace — no silent swallowing,
    # matching the R_optimized fix (commit e4842eb). Swallowing here caused a
    # silent zero-fill of the radar metrics in the negative/permutation controls.
    for name, compute_fn in modules:
        result = compute_fn()
        for sig in names_sigs:
            for ds in names_datasets:
                radar_values[sig][ds].update(result["radar_values"][sig][ds])

    radar_result = compute_radar(radar_values, names_sigs, names_datasets)
    return {
        "radar_values": radar_values,
        "output_table": radar_result["output_table"],
    }


def _run_single_nc_iteration(
    expr_df: pd.DataFrame,
    sig_len: int,
    seed_entropy: int,
    ds_name: str,
    sig_name: str,
) -> np.ndarray:
    """Run one negative control iteration — self-contained for joblib.

    Returns a 1D array of 14 metric values.
    """
    rng = np.random.default_rng(seed_entropy)
    all_genes = list(expr_df.index)
    random_idx = rng.choice(len(all_genes), size=sig_len, replace=False)
    random_genes = [all_genes[i] for i in random_idx]

    random_sig = {sig_name: random_genes}
    result = _compute_qc_metrics(random_sig, [sig_name], {ds_name: expr_df}, [ds_name])

    # Empty output_table indicates a real bug in _compute_qc_metrics — fail loudly.
    if result["output_table"].shape[0] == 0:
        raise RuntimeError(
            f"_compute_qc_metrics returned empty output_table for NC iteration "
            f"(sig={sig_name}, ds={ds_name})"
        )
    return result["output_table"].iloc[0].values


def _run_single_perm_iteration(
    expr_df: pd.DataFrame,
    genes_present: list[str],
    gene_list: list[str],
    seed_entropy: int,
    ds_name: str,
    sig_name: str,
) -> np.ndarray:
    """Run one permutation control iteration — self-contained for joblib.

    Permutes only the signature gene rows, then computes metrics.
    Returns a 1D array of 14 metric values.
    """
    rng = np.random.default_rng(seed_entropy)
    perm_df = expr_df.copy()

    for col in perm_df.columns:
        perm_idx = rng.permutation(len(genes_present))
        perm_df.loc[genes_present, col] = expr_df.loc[
            [genes_present[i] for i in perm_idx], col
        ].values

    perm_sig = {sig_name: gene_list}
    result = _compute_qc_metrics(perm_sig, [sig_name], {ds_name: perm_df}, [ds_name])

    # Empty output_table indicates a real bug in _compute_qc_metrics — fail loudly.
    if result["output_table"].shape[0] == 0:
        raise RuntimeError(
            f"_compute_qc_metrics returned empty output_table for permutation iteration "
            f"(sig={sig_name}, ds={ds_name})"
        )
    return result["output_table"].iloc[0].values


def run_negative_control(
    gene_sigs_list: dict[str, list[str]],
    expression_matrix_list: dict[str, pd.DataFrame],
    num_resampling: int = 50,
    n_jobs: int = -1,
    seed: int = 42,
) -> dict:
    """Run negative and permutation controls with joblib parallelism.

    Args:
        gene_sigs_list: dict of signature name -> gene list
        expression_matrix_list: dict of dataset name -> DataFrame
        num_resampling: number of random resamples
        n_jobs: number of parallel jobs (-1 = all cores, 1 = sequential)
        seed: random seed for reproducibility

    Returns dict with:
        negative_controls: dict [dataset][sig] -> {summary, metrics_table}
        permutation_controls: dict [dataset][sig] -> {summary, metrics_table}
        elapsed_seconds: wall-clock time
    """
    _t0 = time.perf_counter()
    # Generate deterministic, independent seeds for all iterations
    seed_seq = np.random.SeedSequence(seed)

    neg_results: dict = {}
    perm_results: dict = {}

    for ds_name, expr_df in expression_matrix_list.items():
        neg_results[ds_name] = {}
        perm_results[ds_name] = {}

        for sig_name, gene_list in gene_sigs_list.items():
            sig_len = len(gene_list)
            genes_present = [g for g in gene_list if g in expr_df.index]

            # Spawn independent seeds for negative + permutation controls
            nc_seeds = seed_seq.spawn(num_resampling)
            perm_seeds = seed_seq.spawn(num_resampling)

            # --- NEGATIVE CONTROLS (parallel) ---
            nc_metrics_list = Parallel(n_jobs=n_jobs, prefer="processes")(
                delayed(_run_single_nc_iteration)(
                    expr_df, sig_len, nc_seeds[i].entropy, ds_name, f"NC{i+1}"
                )
                for i in range(num_resampling)
            )

            neg_table = pd.DataFrame(
                np.array(nc_metrics_list),
                columns=ALL_METRICS,
                index=[f"NC{i+1}" for i in range(num_resampling)],
            )

            # Vectorized quantiles
            quant_probs = [0.025, 0.25, 0.5, 0.75, 0.975]
            quant_mat = neg_table.quantile(quant_probs)
            mean_vals = neg_table.mean()
            summary = pd.concat([mean_vals.to_frame("mean").T, quant_mat])
            summary.index = ["mean", "Q0.025", "Q0.25", "Q0.5", "Q0.75", "Q0.975"]

            neg_results[ds_name][sig_name] = {
                "summary": summary,
                "metrics_table": neg_table,
            }

            # --- PERMUTATION CONTROLS (parallel) ---
            perm_metrics_list = Parallel(n_jobs=n_jobs, prefer="processes")(
                delayed(_run_single_perm_iteration)(
                    expr_df, genes_present, gene_list,
                    perm_seeds[i].entropy, ds_name, sig_name
                )
                for i in range(num_resampling)
            )

            perm_table = pd.DataFrame(
                np.array(perm_metrics_list),
                columns=ALL_METRICS,
                index=[f"PC{i+1}" for i in range(num_resampling)],
            )

            quant_mat = perm_table.quantile(quant_probs)
            mean_vals = perm_table.mean()
            summary = pd.concat([mean_vals.to_frame("mean").T, quant_mat])
            summary.index = ["mean", "Q0.025", "Q0.25", "Q0.5", "Q0.75", "Q0.975"]

            perm_results[ds_name][sig_name] = {
                "summary": summary,
                "metrics_table": perm_table,
            }

    return {
        "negative_controls": neg_results,
        "permutation_controls": perm_results,
        "elapsed_seconds": time.perf_counter() - _t0,
    }
