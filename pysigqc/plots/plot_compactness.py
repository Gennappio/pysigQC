"""Autocorrelation heatmaps and density overlay plot.

Produces:
- ``sig_autocor_hmps.pdf`` — clustered heatmaps of gene-gene Spearman
  correlations, one per (sig, dataset). Uses seaborn clustermap with
  blue-white-red diverging palette and dendrograms (matching R heatmap.2).
- ``sig_autocor_dens.pdf`` — overlaid density curves of pairwise correlations,
  one line per (sig, dataset). Color=dataset, linestyle=signature.
  Legend positioned outside right (matching R's par(omi=...)).

Matches R's eval_compactness_loc.R styling.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import gaussian_kde

from ._colors import dataset_colors, sig_linestyle
from ._style import save_pdf, dynamic_fontsize


# R's gplots::colorpanel(100,"blue","white","red")
_HEATMAP_CMAP = "RdBu_r"


def plot_compactness(
    compact_result: dict,
    names_sigs: list[str],
    names_datasets: list[str],
    out_dir: str | Path,
) -> list[Path]:
    autocor_matrices = compact_result["autocor_matrices"]
    n_sigs = len(names_sigs)
    n_ds = len(names_datasets)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    max_label = max(len(f"{ds} {sig}") for sig in names_sigs for ds in names_datasets)
    ds_colors = dataset_colors(n_ds)

    paths: list[Path] = []

    # --- 1. Autocorrelation heatmaps (clustermap per sig x dataset) ---
    # Use PdfPages to put all heatmaps in one PDF, matching R's single PDF
    pdf_path = out_dir / "sig_autocor_hmps.pdf"
    with PdfPages(pdf_path) as pp:
        for ki, sig in enumerate(names_sigs):
            for di, ds in enumerate(names_datasets):
                mat = autocor_matrices[sig][ds]
                if mat.shape[0] <= 1:
                    continue
                n_genes = mat.shape[0]
                # Adaptive label size: R's max(min(0.5, 4*4/n), 0.06)
                label_fs = max(min(5, 4 * 4 / n_genes * 10), 3)
                # Adaptive margins: R's 1 + max(nchar(names))/2
                max_name_len = max(len(str(g)) for g in mat.index) if n_genes > 0 else 5
                cell_size = max(0.3, min(0.6, 4.0 / n_genes))
                fig_size = max(4, n_genes * cell_size + 2)

                try:
                    g = sns.clustermap(
                        mat.astype(float),
                        vmin=-1, vmax=1,
                        cmap=_HEATMAP_CMAP,
                        center=0,
                        row_cluster=True, col_cluster=True,
                        dendrogram_ratio=(0.15, 0.15),
                        figsize=(fig_size, fig_size),
                        xticklabels=True, yticklabels=True,
                        linewidths=0.3, linecolor="white",
                        cbar_kws={"label": "Rho", "shrink": 0.6},
                    )
                    g.ax_heatmap.tick_params(labelsize=label_fs)
                    g.ax_heatmap.set_title(f"{ds}  {sig}", fontsize=9, pad=10)
                    # Grey for NA values
                    g.ax_heatmap.set_facecolor("#D3D3D3")
                    pp.savefig(g.fig, bbox_inches="tight")
                    plt.close(g.fig)
                except Exception:
                    pass
    paths.append(pdf_path)

    # --- 2. Autocorrelation density overlay ---
    # R: single 5x5 inch plot, legend outside right, par(omi=c(0,0,0,w))
    fig, ax = plt.subplots(figsize=(5, 5))
    legend_handles = []
    legend_labels = []

    # First pass: compute max density for ylim (matching R's approach)
    max_dens = 0
    for ki, sig in enumerate(names_sigs):
        for di, ds in enumerate(names_datasets):
            mat = autocor_matrices[sig][ds]
            if mat.shape[0] <= 1:
                continue
            vals = mat.values[np.triu_indices_from(mat.values, k=1)]
            vals = vals[np.isfinite(vals)]
            if len(vals) >= 2:
                try:
                    dens_max = gaussian_kde(vals)(np.linspace(-1, 1, 300)).max()
                    max_dens = max(max_dens, dens_max)
                except np.linalg.LinAlgError:
                    pass

    # Second pass: draw
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

    ax.set_xlim(-1, 1)
    if max_dens > 0:
        ax.set_ylim(0, np.ceil(max_dens))
    ax.set_xlabel("Rho", fontsize=10)
    ax.set_ylabel("Density", fontsize=10)
    ax.set_title("Intra-sig. Corr. Density", fontsize=11)

    # Legend outside right (R: legend at par('usr')[2]+0.05, par('usr')[4])
    if legend_handles:
        legend_font = min(5.0, 4 * 10 / max_label)
        ax.legend(legend_handles, legend_labels, fontsize=legend_font,
                  loc="center left", bbox_to_anchor=(1.02, 0.5),
                  borderaxespad=0, frameon=False)

    fig.tight_layout()
    paths.append(save_pdf(fig, out_dir, "sig_autocor_dens.pdf"))

    return paths
