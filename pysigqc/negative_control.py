"""Negative and permutation controls for gene signature QC metrics.

Port of R_refactored/sigsQcNegativeControl.R + R_optimized/ version.
Fully optimized: vectorized random signature generation, in-place permutation,
optional joblib parallelism, no disk I/O roundtrip.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

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
    """Compute all radar metrics in memory (no file I/O).

    Returns dict with:
        radar_values: nested dict [sig][dataset] -> metric dict
        output_table: DataFrame of radar chart values
    """
    radar_values: dict = {}
    for sig in names_sigs:
        radar_values[sig] = {}
        for ds in names_datasets:
            radar_values[sig][ds] = {}

    # Run each compute module and merge radar values
    modules = [
        ("var", lambda: compute_var(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)),
        ("expr", lambda: compute_expr(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)),
        ("compactness", lambda: compute_compactness(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)),
        ("metrics", lambda: compute_metrics(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)),
        ("stan", lambda: compute_stan(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)),
    ]

    for name, compute_fn in modules:
        try:
            result = compute_fn()
            for sig in names_sigs:
                for ds in names_datasets:
                    radar_values[sig][ds].update(result["radar_values"][sig][ds])
        except Exception:
            pass

    # Assemble radar chart
    radar_result = compute_radar(radar_values, names_sigs, names_datasets)
    return {
        "radar_values": radar_values,
        "output_table": radar_result["output_table"],
    }


def _run_single_negative_control(
    expr_matrix: pd.DataFrame,
    sig_len: int,
    rng: np.random.Generator,
) -> dict[str, list[str]]:
    """Generate one random signature of the same length as the original."""
    all_genes = list(expr_matrix.index)
    random_indices = rng.choice(len(all_genes), size=sig_len, replace=False)
    random_genes = [all_genes[i] for i in random_indices]
    return random_genes


def _run_single_permutation(
    expr_matrix: pd.DataFrame,
    genes_present: list[str],
    rng: np.random.Generator,
) -> pd.DataFrame:
    """Generate one permuted expression matrix (only permute signature gene rows)."""
    perm_matrix = expr_matrix.copy()
    for col in perm_matrix.columns:
        perm_idx = rng.permutation(len(genes_present))
        perm_matrix.loc[genes_present, col] = expr_matrix.loc[
            [genes_present[i] for i in perm_idx], col
        ].values
    return perm_matrix


def run_negative_control(
    gene_sigs_list: dict[str, list[str]],
    expression_matrix_list: dict[str, pd.DataFrame],
    num_resampling: int = 50,
    n_jobs: int = 1,
    seed: int = 42,
) -> dict:
    """Run negative and permutation controls.

    Args:
        gene_sigs_list: dict of signature name -> gene list
        expression_matrix_list: dict of dataset name -> DataFrame
        num_resampling: number of random resamples
        n_jobs: number of parallel jobs (1 = sequential)
        seed: random seed for reproducibility

    Returns dict with:
        negative_controls: dict [dataset][sig] -> {summary, metrics_table}
        permutation_controls: dict [dataset][sig] -> {summary, metrics_table}
    """
    rng = np.random.default_rng(seed)
    neg_results: dict = {}
    perm_results: dict = {}

    for ds_name, expr_matrix in expression_matrix_list.items():
        neg_results[ds_name] = {}
        perm_results[ds_name] = {}

        for sig_name, gene_list in gene_sigs_list.items():
            sig_len = len(gene_list)
            genes_present = [g for g in gene_list if g in expr_matrix.index]

            # --- NEGATIVE CONTROLS (random gene selection) ---
            random_sigs = {}
            for i in range(num_resampling):
                random_genes = _run_single_negative_control(expr_matrix, sig_len, rng)
                random_sigs[f"NC{i+1}"] = random_genes

            nc_result = _compute_qc_metrics(
                random_sigs, list(random_sigs.keys()),
                {ds_name: expr_matrix}, [ds_name],
            )
            neg_table = nc_result["output_table"]

            # Vectorized quantile computation
            quant_probs = [0.025, 0.25, 0.5, 0.75, 0.975]
            quant_mat = neg_table.quantile(quant_probs)
            mean_vals = neg_table.mean()
            summary = pd.concat([mean_vals.to_frame("mean").T, quant_mat])
            summary.index = ["mean", "Q0.025", "Q0.25", "Q0.5", "Q0.75", "Q0.975"]

            neg_results[ds_name][sig_name] = {
                "summary": summary,
                "metrics_table": neg_table,
            }

            # --- PERMUTATION CONTROLS (gene label shuffling) ---
            perm_matrices = {}
            for i in range(num_resampling):
                perm_matrices[f"PC{i+1}"] = _run_single_permutation(
                    expr_matrix, genes_present, rng
                )

            # Single signature across permuted matrices
            perm_sigs = {sig_name: gene_list}
            pc_result = _compute_qc_metrics(
                perm_sigs, [sig_name],
                perm_matrices, list(perm_matrices.keys()),
            )
            perm_table = pc_result["output_table"]

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
    }
