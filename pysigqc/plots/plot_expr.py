"""Expression level plots matching R's eval_expr_loc.R styling.

Produces:
- ``sig_expr_barcharts_NA_values.pdf`` — sorted descending NA proportions
- ``sig_expr_barcharts.pdf`` — sorted ascending expression proportions
- ``sig_expr_density_plots.pdf`` — KDE of expression proportions (adjust=0.25)

R style: barplot with border=NA, gene names at 45 degrees via text(),
adaptive font sizing, y-axis drawn separately.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ._style import figure_grid, save_pdf, dynamic_fontsize, gene_label_fontsize


def _draw_barchart(ax, vals, title, ylabel, font, ascending=False):
    """Draw a single bar chart matching R's barplot() style."""
    if ascending:
        vals = vals.sort_values(ascending=True)
    else:
        vals = vals.sort_values(ascending=False)

    n = len(vals)
    label_fs = gene_label_fontsize(n)

    bars = ax.bar(range(n), vals.values, width=0.8,
                  color="#D3D3D3", edgecolor="none")  # R default grey, border=NA

    # Gene names at 45 degrees (R: text(srt=45, adj=c(1.1,1.1)))
    ax.set_xticks(range(n))
    ax.set_xticklabels(vals.index, rotation=45, ha="right",
                       fontsize=label_fs)

    ax.set_ylim(0, 1)
    ax.set_title(title, fontsize=font)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.tick_params(axis="y", labelsize=6)


def plot_expr(
    expr_result: dict,
    names_sigs: list[str],
    names_datasets: list[str],
    out_dir: str | Path,
) -> list[Path]:
    na_props = expr_result["na_proportions"]
    expr_props = expr_result["expr_proportions"]

    n_sigs = len(names_sigs)
    n_ds = len(names_datasets)

    max_label = max(len(f"{ds} {sig}") for sig in names_sigs for ds in names_datasets)
    font = min(10.0, 4 * 12 / max_label)  # R: min(1, 4*12/max_line_length) * 10

    paths: list[Path] = []

    # --- 1. NA proportion barcharts (sorted descending) ---
    fig, axes = figure_grid(n_sigs, n_ds)
    for ki, sig in enumerate(names_sigs):
        for di, ds in enumerate(names_datasets):
            _draw_barchart(
                axes[ki, di], na_props[sig][ds],
                f"Proportion NA values\n{ds} {sig}",
                "Proportion NA", font, ascending=False,
            )
    fig.tight_layout()
    paths.append(save_pdf(fig, out_dir, "sig_expr_barcharts_NA_values.pdf"))

    # --- 2. Expression proportion barcharts (sorted ascending) ---
    fig, axes = figure_grid(n_sigs, n_ds)
    for ki, sig in enumerate(names_sigs):
        for di, ds in enumerate(names_datasets):
            _draw_barchart(
                axes[ki, di], expr_props[sig][ds],
                f"Proportion expressed\n{ds} {sig}",
                "Proportion above threshold", font, ascending=True,
            )
    fig.tight_layout()
    paths.append(save_pdf(fig, out_dir, "sig_expr_barcharts.pdf"))

    # --- 3. Expression density plots (R: density(adjust=0.25)) ---
    fig, axes = figure_grid(n_sigs, n_ds)
    for ki, sig in enumerate(names_sigs):
        for di, ds in enumerate(names_datasets):
            ax = axes[ki, di]
            vals = expr_props[sig][ds].dropna().values
            if len(vals) > 1:
                from scipy.stats import gaussian_kde
                try:
                    # R adjust=0.25 ≈ bw_method * 0.25
                    kde = gaussian_kde(vals, bw_method=0.25)
                    x = np.linspace(0, 1, 200)
                    y = kde(x)
                    ax.plot(x, y, color="black", lw=1.5)
                    ax.fill_between(x, y, alpha=0.1, color="grey")
                except np.linalg.LinAlgError:
                    ax.hist(vals, bins=20, density=True, color="#D3D3D3",
                            edgecolor="none")
            ax.set_xlim(0, 1)
            ax.set_title(f"Density of expr. proportion\n{ds} {sig}", fontsize=font)
            ax.set_xlabel("Proportion above threshold", fontsize=8)
            ax.set_ylabel("Density", fontsize=8)
    fig.tight_layout()
    paths.append(save_pdf(fig, out_dir, "sig_expr_density_plots.pdf"))

    return paths
