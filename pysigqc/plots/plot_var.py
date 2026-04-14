"""Mean vs SD scatter plot with quantile reference lines.

Produces ``sig_mean_vs_sd.pdf`` — one subplot per (signature, dataset) pair.
Grey dots = all genes, red dots = signature genes, dashed lines at quantiles.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ._colors import BACKGROUND_COLOR, SIG_GENE_COLOR, QUANTILE_COLORS
from ._style import figure_grid, save_pdf, dynamic_fontsize


def plot_var(
    var_result: dict,
    names_sigs: list[str],
    names_datasets: list[str],
    out_dir: str | Path,
) -> Path:
    """Generate Mean vs SD scatter. Returns path to the PDF."""
    mean_sd_tables = var_result["mean_sd_tables"]
    all_sd = var_result["all_sd"]
    all_mean = var_result["all_mean"]
    inter = var_result["inter"]

    n_sigs = len(names_sigs)
    n_ds = len(names_datasets)

    max_label = max(len(f"{ds} {sig}") for sig in names_sigs for ds in names_datasets)
    font = dynamic_fontsize(max_label)

    fig, axes = figure_grid(n_sigs, n_ds)

    for ki, sig in enumerate(names_sigs):
        for di, ds in enumerate(names_datasets):
            ax = axes[ki, di]

            sd_all = all_sd[sig][ds].dropna().values
            mean_all = all_mean[sig][ds].dropna().values
            genes = inter[sig][ds]

            # All genes — grey
            ax.scatter(mean_all, sd_all, s=3, c=BACKGROUND_COLOR, alpha=0.6, edgecolors="none")

            # Signature genes — red
            tbl = mean_sd_tables[sig][ds]
            ax.scatter(tbl["Mean"].values, tbl["SD"].values, s=6, c=SIG_GENE_COLOR, edgecolors="none", zorder=3)

            # Quantile lines for mean (vertical) and SD (horizontal)
            finite_mean = mean_all[np.isfinite(mean_all)]
            finite_sd = sd_all[np.isfinite(sd_all)]
            if len(finite_mean) > 0 and len(finite_sd) > 0:
                qm = np.quantile(finite_mean, [0.1, 0.25, 0.5, 0.75, 0.9])
                qs = np.quantile(finite_sd, [0.1, 0.25, 0.5, 0.75, 0.9])
                for idx, color_key in zip([0, 4], ["outer"] * 2):
                    ax.axvline(qm[idx], ls="--", lw=0.7, color=QUANTILE_COLORS["outer"])
                    ax.axhline(qs[idx], ls="--", lw=0.7, color=QUANTILE_COLORS["outer"])
                for idx in [1, 3]:
                    ax.axvline(qm[idx], ls="--", lw=0.7, color=QUANTILE_COLORS["inner"])
                    ax.axhline(qs[idx], ls="--", lw=0.7, color=QUANTILE_COLORS["inner"])
                ax.axvline(qm[2], ls="--", lw=0.7, color=QUANTILE_COLORS["median"])
                ax.axhline(qs[2], ls="--", lw=0.7, color=QUANTILE_COLORS["median"])

            ax.set_title(f"Mean vs SD\n{ds} {sig}", fontsize=font)
            ax.set_xlabel("Mean", fontsize=8)
            ax.set_ylabel("Standard deviation", fontsize=8)

    fig.tight_layout()
    return save_pdf(fig, out_dir, "sig_mean_vs_sd.pdf")
