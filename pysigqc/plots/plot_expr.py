"""Expression level plots: NA proportion barcharts, expression proportion
barcharts, and expression density plots.

Produces three PDFs matching R's eval_expr_loc output.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ._style import figure_grid, save_pdf, dynamic_fontsize


def plot_expr(
    expr_result: dict,
    names_sigs: list[str],
    names_datasets: list[str],
    out_dir: str | Path,
) -> list[Path]:
    """Generate expression level plots. Returns list of PDF paths."""
    na_props = expr_result["na_proportions"]
    expr_props = expr_result["expr_proportions"]

    n_sigs = len(names_sigs)
    n_ds = len(names_datasets)

    max_label = max(len(f"{ds} {sig}") for sig in names_sigs for ds in names_datasets)
    font = dynamic_fontsize(max_label)

    paths: list[Path] = []

    # --- 1. NA proportion barcharts ---
    fig, axes = figure_grid(n_sigs, n_ds)
    for ki, sig in enumerate(names_sigs):
        for di, ds in enumerate(names_datasets):
            ax = axes[ki, di]
            vals = na_props[sig][ds].sort_values(ascending=False)
            ax.bar(range(len(vals)), vals.values, color="#808080", edgecolor="none")
            ax.set_xticks(range(len(vals)))
            ax.set_xticklabels(vals.index, rotation=90, fontsize=max(4, font * 0.6))
            ax.set_ylim(0, 1)
            ax.set_title(f"Prop. NA\n{ds} {sig}", fontsize=font)
            ax.set_ylabel("Proportion NA", fontsize=8)
    fig.tight_layout()
    paths.append(save_pdf(fig, out_dir, "sig_expr_barcharts_NA_values.pdf"))

    # --- 2. Expression proportion barcharts ---
    fig, axes = figure_grid(n_sigs, n_ds)
    for ki, sig in enumerate(names_sigs):
        for di, ds in enumerate(names_datasets):
            ax = axes[ki, di]
            vals = expr_props[sig][ds].sort_values(ascending=True)
            ax.bar(range(len(vals)), vals.values, color="#808080", edgecolor="none")
            ax.set_xticks(range(len(vals)))
            ax.set_xticklabels(vals.index, rotation=90, fontsize=max(4, font * 0.6))
            ax.set_ylim(0, 1)
            ax.set_title(f"Prop. expressed\n{ds} {sig}", fontsize=font)
            ax.set_ylabel("Proportion above threshold", fontsize=8)
    fig.tight_layout()
    paths.append(save_pdf(fig, out_dir, "sig_expr_barcharts.pdf"))

    # --- 3. Expression density plots ---
    fig, axes = figure_grid(n_sigs, n_ds)
    for ki, sig in enumerate(names_sigs):
        for di, ds in enumerate(names_datasets):
            ax = axes[ki, di]
            vals = expr_props[sig][ds].dropna().values
            if len(vals) > 1:
                from scipy.stats import gaussian_kde
                try:
                    kde = gaussian_kde(vals, bw_method=0.25)
                    x = np.linspace(0, 1, 200)
                    ax.plot(x, kde(x), color="black", lw=1.5)
                    ax.fill_between(x, kde(x), alpha=0.15, color="black")
                except np.linalg.LinAlgError:
                    ax.hist(vals, bins=20, density=True, color="#808080")
            ax.set_xlim(0, 1)
            ax.set_title(f"Expression density\n{ds} {sig}", fontsize=font)
            ax.set_xlabel("Proportion above threshold", fontsize=8)
            ax.set_ylabel("Density", fontsize=8)
    fig.tight_layout()
    paths.append(save_pdf(fig, out_dir, "sig_expr_density_plots.pdf"))

    return paths
