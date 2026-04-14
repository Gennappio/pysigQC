"""Standardization comparison scatter plot.

Produces ``sig_standardisation_comp.pdf`` — one subplot per (sig, dataset).
Each subplot is a density-colored scatter of raw median score vs z-transformed
median score, with a marginal density overlay in red.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ._colors import JET_CMAP
from ._style import figure_grid, save_pdf, dynamic_fontsize


def plot_stan(
    stan_result: dict,
    names_sigs: list[str],
    names_datasets: list[str],
    out_dir: str | Path,
) -> Path:
    """Generate standardization comparison scatter. Returns PDF path."""
    med_scores = stan_result["med_scores"]
    z_scores = stan_result["z_transf_scores"]

    n_sigs = len(names_sigs)
    n_ds = len(names_datasets)

    max_label = max(len(f"{ds} {sig}") for sig in names_sigs for ds in names_datasets)
    font = dynamic_fontsize(max_label)

    fig, axes = figure_grid(n_sigs, n_ds)

    for ki, sig in enumerate(names_sigs):
        for di, ds in enumerate(names_datasets):
            ax = axes[ki, di]
            x = np.asarray(med_scores[sig][ds], dtype=float)
            y = np.asarray(z_scores[sig][ds], dtype=float)

            mask = np.isfinite(x) & np.isfinite(y)
            x, y = x[mask], y[mask]

            if len(x) > 2:
                ax.hexbin(x, y, gridsize=30, cmap=JET_CMAP, mincnt=1)

                # Marginal density overlay
                from scipy.stats import gaussian_kde, spearmanr
                try:
                    kde = gaussian_kde(x)
                    xs = np.linspace(x.min(), x.max(), 200)
                    dens = kde(xs)
                    # Scale density into upper 30% of y-axis
                    ymin, ymax = ax.get_ylim()
                    dens_scaled = ymin + (ymax - ymin) * 0.7 + dens / dens.max() * (ymax - ymin) * 0.25
                    ax.plot(xs, dens_scaled, color="red", lw=1)
                except np.linalg.LinAlgError:
                    pass

                rho, _ = spearmanr(x, y)
                ax.text(0.05, 0.95, f"rho={rho:.3f}", transform=ax.transAxes,
                        fontsize=7, va="top", ha="left",
                        bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7))
            else:
                ax.text(0.5, 0.5, "Insufficient data", transform=ax.transAxes,
                        ha="center", va="center", fontsize=8)

            ax.set_title(f"Standardisation\n{ds} {sig}", fontsize=font)
            ax.set_xlabel("Median score (raw)", fontsize=8)
            ax.set_ylabel("Z-median score", fontsize=8)

    fig.tight_layout()
    return save_pdf(fig, out_dir, "sig_standardisation_comp.pdf")
