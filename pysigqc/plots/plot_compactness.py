"""Autocorrelation heatmaps and density overlay plot.

Produces:
- ``sig_autocor_hmps.pdf`` — clustered heatmap of gene-gene Spearman
  correlations per (sig, dataset)
- ``sig_autocor_dens.pdf`` — overlaid density curves of pairwise correlations,
  one line per (sig, dataset), colored by dataset, styled by signature
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from ._colors import HEATMAP_CMAP, dataset_colors, sig_linestyle
from ._style import figure_grid, save_pdf, dynamic_fontsize


def plot_compactness(
    compact_result: dict,
    names_sigs: list[str],
    names_datasets: list[str],
    out_dir: str | Path,
) -> list[Path]:
    """Generate compactness plots. Returns list of PDF paths."""
    autocor_matrices = compact_result["autocor_matrices"]
    n_sigs = len(names_sigs)
    n_ds = len(names_datasets)

    max_label = max(len(f"{ds} {sig}") for sig in names_sigs for ds in names_datasets)
    font = dynamic_fontsize(max_label)
    ds_colors = dataset_colors(n_ds)

    paths: list[Path] = []

    # --- 1. Autocorrelation heatmaps (one subplot per sig x dataset) ---
    fig, axes = figure_grid(n_sigs, n_ds, cell_w=5, cell_h=5)
    for ki, sig in enumerate(names_sigs):
        for di, ds in enumerate(names_datasets):
            ax = axes[ki, di]
            mat = autocor_matrices[sig][ds]
            if mat.shape[0] > 1:
                sns.heatmap(
                    mat, ax=ax, vmin=-1, vmax=1, cmap=HEATMAP_CMAP,
                    square=True, linewidths=0.5,
                    xticklabels=True, yticklabels=True,
                    cbar_kws={"shrink": 0.6},
                )
                ax.tick_params(labelsize=max(4, font * 0.5))
            else:
                ax.text(0.5, 0.5, "< 2 genes", transform=ax.transAxes,
                        ha="center", va="center")
            ax.set_title(f"Autocorrelation\n{ds} {sig}", fontsize=font)
    fig.tight_layout()
    paths.append(save_pdf(fig, out_dir, "sig_autocor_hmps.pdf"))

    # --- 2. Autocorrelation density overlay ---
    fig, ax = plt.subplots(figsize=(6, 5))
    from scipy.stats import gaussian_kde
    legend_handles = []
    legend_labels = []

    for ki, sig in enumerate(names_sigs):
        for di, ds in enumerate(names_datasets):
            mat = autocor_matrices[sig][ds]
            if mat.shape[0] <= 1:
                continue
            vals = mat.values[np.triu_indices_from(mat.values, k=1)]
            vals = vals[np.isfinite(vals)]
            if len(vals) < 2:
                continue
            try:
                kde = gaussian_kde(vals)
                x = np.linspace(-1, 1, 300)
                line, = ax.plot(x, kde(x), color=ds_colors[di],
                                linestyle=sig_linestyle(ki), lw=2)
                legend_handles.append(line)
                legend_labels.append(f"{ds} {sig}")
            except np.linalg.LinAlgError:
                pass

    ax.set_xlabel("Rho", fontsize=10)
    ax.set_ylabel("Density", fontsize=10)
    ax.set_title("Intra-sig. Corr. Density", fontsize=12)
    ax.set_xlim(-1, 1)
    if legend_handles:
        ax.legend(legend_handles, legend_labels, fontsize=7,
                  loc="upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0)
    fig.tight_layout()
    paths.append(save_pdf(fig, out_dir, "sig_autocor_dens.pdf"))

    return paths
