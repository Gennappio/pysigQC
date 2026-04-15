"""Mean vs SD scatter plot with quantile reference lines.

Produces ``sig_mean_vs_sd.pdf`` — grid of subplots (rows=sigs, cols=datasets).
Grey filled circles = all genes, red filled circles = signature genes.
Dashed quantile lines at 10/25/50/75/90 percentiles for both axes.
Matches R's eval_var_loc.R styling.
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
    mean_sd_tables = var_result["mean_sd_tables"]
    all_sd = var_result["all_sd"]
    all_mean = var_result["all_mean"]
    inter = var_result["inter"]

    n_sigs = len(names_sigs)
    n_ds = len(names_datasets)

    max_label = max(len(f"{ds} {sig}") for sig in names_sigs for ds in names_datasets)
    title_cex = min(1.0, 4 * 10 / max_label)  # R's cex.main formula
    font = title_cex * 10  # scale to matplotlib points

    fig, axes = figure_grid(n_sigs, n_ds)

    for ki, sig in enumerate(names_sigs):
        for di, ds in enumerate(names_datasets):
            ax = axes[ki, di]

            sd_all = all_sd[sig][ds].dropna().values
            mean_all = all_mean[sig][ds].dropna().values

            # All genes — grey filled circles (R: pch=19, col='grey', cex=0.5)
            ax.scatter(mean_all, sd_all, s=8, c="grey", edgecolors="none",
                       alpha=0.7, zorder=1, rasterized=True)

            # Signature genes — red filled circles (R: pch=19, col='red', cex=0.5)
            tbl = mean_sd_tables[sig][ds]
            ax.scatter(tbl["Mean"].values, tbl["SD"].values,
                       s=8, c="red", edgecolors="none", zorder=3)

            # Quantile lines (R: lty=5 = long dash)
            finite_mean = mean_all[np.isfinite(mean_all)]
            finite_sd = sd_all[np.isfinite(sd_all)]
            if len(finite_mean) > 0 and len(finite_sd) > 0:
                qm = np.quantile(finite_mean, [0.1, 0.25, 0.5, 0.75, 0.9])
                qs = np.quantile(finite_sd, [0.1, 0.25, 0.5, 0.75, 0.9])
                dash = (8, 4)  # R lty=5 (long dash)
                for idx in (0, 4):  # 10%, 90%
                    ax.axvline(qm[idx], ls="--", dashes=dash, lw=0.8,
                               color=QUANTILE_COLORS["outer"])
                    ax.axhline(qs[idx], ls="--", dashes=dash, lw=0.8,
                               color=QUANTILE_COLORS["outer"])
                for idx in (1, 3):  # 25%, 75%
                    ax.axvline(qm[idx], ls="--", dashes=dash, lw=0.8,
                               color=QUANTILE_COLORS["inner"])
                    ax.axhline(qs[idx], ls="--", dashes=dash, lw=0.8,
                               color=QUANTILE_COLORS["inner"])
                # 50% (median)
                ax.axvline(qm[2], ls="--", dashes=dash, lw=0.8,
                           color=QUANTILE_COLORS["median"])
                ax.axhline(qs[2], ls="--", dashes=dash, lw=0.8,
                           color=QUANTILE_COLORS["median"])

            # Legend matching R: topright, no box, cex=0.5
            from matplotlib.lines import Line2D
            handles = [
                Line2D([0], [0], color=QUANTILE_COLORS["outer"], ls="--", lw=0.8),
                Line2D([0], [0], color=QUANTILE_COLORS["inner"], ls="--", lw=0.8),
                Line2D([0], [0], color=QUANTILE_COLORS["median"], ls="--", lw=0.8),
            ]
            ax.legend(handles, ["10%, 90%", "25%, 75%", "50%"],
                      loc="upper right", fontsize=5, frameon=False,
                      handlelength=1.5)

            ax.set_title(f"Mean vs SD of expression\n{ds} {sig}", fontsize=font)
            ax.set_xlabel("Mean", fontsize=8)
            ax.set_ylabel("Standard deviation", fontsize=8)

    fig.tight_layout()
    return save_pdf(fig, out_dir, "sig_mean_vs_sd.pdf")
