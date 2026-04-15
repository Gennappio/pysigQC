"""Standardization comparison: raw median score vs z-transformed median score.

Produces ``sig_standardisation_comp.pdf`` — grid (rows=sigs, cols=datasets).
Each subplot: smoothScatter-style 2D KDE scatter with jet colormap,
red density overlay on a secondary y-axis (right), Spearman rho annotation.
Matches R's eval_stan_loc.R styling.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import gaussian_kde, spearmanr

from ._style import figure_grid, save_pdf, dynamic_fontsize

# R's jet.colors: colorRampPalette(c("#00007F","blue","#007FFF","cyan",
#                                     "#7FFF7F","yellow","#FF7F00","red","#7F0000"))
_JET_COLORS = ["#00007F", "#0000FF", "#007FFF", "#00FFFF",
               "#7FFF7F", "#FFFF00", "#FF7F00", "#FF0000", "#7F0000"]
_JET_CMAP = LinearSegmentedColormap.from_list("r_jet", _JET_COLORS, N=256)


def _smoothscatter(ax, x, y, title, xlabel, ylabel, font):
    """R-style smoothScatter: 2D KDE coloured scatter + red density overlay."""
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]

    if len(x) < 3:
        ax.text(0.5, 0.5, "Insufficient data", transform=ax.transAxes,
                ha="center", va="center", fontsize=8, color="grey")
        ax.set_title(title, fontsize=font)
        return

    # 2D KDE density-coloured scatter (matches smoothScatter)
    try:
        xy = np.vstack([x, y])
        kde = gaussian_kde(xy)
        density = kde(xy)
        idx = density.argsort()  # plot low-density first
        ax.scatter(x[idx], y[idx], c=density[idx], s=4, cmap=_JET_CMAP,
                   edgecolors="none", rasterized=True)
    except np.linalg.LinAlgError:
        ax.scatter(x, y, s=4, c="blue", edgecolors="none", rasterized=True)

    # Red density overlay on right y-axis (R: lines(density(x), col='red', lwd=2))
    try:
        kde_1d = gaussian_kde(x)
        xs = np.linspace(x.min(), x.max(), 200)
        dens = kde_1d(xs)
        ax2 = ax.twinx()
        ax2.plot(xs, dens, color="red", lw=2, alpha=0.8)
        ax2.set_ylabel("Density", fontsize=7, color="red")
        ax2.tick_params(axis="y", labelsize=5, colors="red")
    except np.linalg.LinAlgError:
        pass

    # Spearman rho annotation
    rho, _ = spearmanr(x, y)
    ax.set_title(f"{title}\nrho = {rho:.3f}", fontsize=font)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)


def plot_stan(
    stan_result: dict,
    names_sigs: list[str],
    names_datasets: list[str],
    out_dir: str | Path,
) -> Path:
    med_scores = stan_result["med_scores"]
    z_scores = stan_result["z_transf_scores"]

    n_sigs = len(names_sigs)
    n_ds = len(names_datasets)

    max_label = max(len(f"{ds} {sig}") for sig in names_sigs for ds in names_datasets)
    font = min(10.0, 4 * 10 / max_label)

    fig, axes = figure_grid(n_sigs, n_ds)

    for ki, sig in enumerate(names_sigs):
        for di, ds in enumerate(names_datasets):
            x = np.asarray(med_scores[sig][ds], dtype=float)
            y = np.asarray(z_scores[sig][ds], dtype=float)
            _smoothscatter(axes[ki, di], x, y,
                           f"Median vs Z-median\n{ds} {sig}",
                           "Median score", "Z-median score", font)

    fig.tight_layout()
    return save_pdf(fig, out_dir, "sig_standardisation_comp.pdf")
