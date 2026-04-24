"""Autocorrelation heatmaps and density overlay plot.

Produces:
- ``sig_autocor_hmps.pdf`` — clustered heatmaps of gene-gene Spearman
  correlations, one per (sig, dataset). Uses seaborn clustermap with
  blue-white-red diverging palette and dendrograms (matching R heatmap.2).
- ``sig_autocor_dens.pdf`` — overlaid density curves of pairwise correlations,
  one line per (sig, dataset). Color=dataset, linestyle=signature.
  Legend positioned outside right (matching R's par(omi=...)).
- ``sig_autocor_rankProd_{sig}.pdf`` — RankProd plots showing genes with
  consistently high/low autocorrelation across datasets (matching R's plotRP).

Matches R's eval_compactness_loc.R styling.
"""

from __future__ import annotations

import logging
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

logger = logging.getLogger(__name__)


# R's gplots::colorpanel(100,"blue","white","red")
_HEATMAP_CMAP = "RdBu_r"


def _plot_rank_product(
    rp_table,
    sig_name: str,
    out_dir: Path,
) -> Path | None:
    """Generate RankProd plot matching R's RankProd::plotRP().

    Creates a 1x2 subplot:
      - Left: pfp_up vs gene rank (genes with consistently HIGH autocorrelation)
      - Right: pfp_down vs gene rank (genes with consistently LOW autocorrelation)

    Args:
        rp_table: DataFrame with columns pfp_up, pfp_down, rp_up, rp_down
        sig_name: signature name for title and filename
        out_dir: output directory

    Returns:
        Path to generated PDF, or None if no data
    """
    import pandas as pd

    if rp_table is None or (isinstance(rp_table, pd.DataFrame) and rp_table.empty):
        return None

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    n_genes = len(rp_table)

    # --- Left plot: Up-regulated (high autocorrelation) ---
    ax = axes[0]
    sorted_up = rp_table.sort_values("rp_up")
    ranks = np.arange(1, n_genes + 1)
    pfp_up = sorted_up["pfp_up"].values

    # Plot points
    ax.scatter(pfp_up, ranks, c="blue", s=20, alpha=0.7, edgecolors="none")

    # Add gene labels for top genes (pfp < 0.05)
    sig_genes_up = sorted_up[sorted_up["pfp_up"] < 0.05]
    for gene in sig_genes_up.index:
        rank_pos = np.where(sorted_up.index == gene)[0][0] + 1
        pfp_val = sig_genes_up.loc[gene, "pfp_up"]
        ax.annotate(
            str(gene), (pfp_val, rank_pos),
            fontsize=6, ha="left", va="center",
            xytext=(3, 0), textcoords="offset points"
        )

    # Add vertical line at pfp=0.05
    ax.axvline(x=0.05, color="red", linestyle="--", linewidth=1, alpha=0.7)

    ax.set_xlabel("Estimated PFP", fontsize=10)
    ax.set_ylabel("Number of identified genes", fontsize=10)
    ax.set_title(f"Class 1: High autocorrelation\n{sig_name}", fontsize=10)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(0, n_genes + 1)

    # --- Right plot: Down-regulated (low autocorrelation) ---
    ax = axes[1]
    sorted_down = rp_table.sort_values("rp_down")
    pfp_down = sorted_down["pfp_down"].values

    ax.scatter(pfp_down, ranks, c="blue", s=20, alpha=0.7, edgecolors="none")

    # Add gene labels for top genes (pfp < 0.05)
    sig_genes_down = sorted_down[sorted_down["pfp_down"] < 0.05]
    for gene in sig_genes_down.index:
        rank_pos = np.where(sorted_down.index == gene)[0][0] + 1
        pfp_val = sig_genes_down.loc[gene, "pfp_down"]
        ax.annotate(
            str(gene), (pfp_val, rank_pos),
            fontsize=6, ha="left", va="center",
            xytext=(3, 0), textcoords="offset points"
        )

    ax.axvline(x=0.05, color="red", linestyle="--", linewidth=1, alpha=0.7)

    ax.set_xlabel("Estimated PFP", fontsize=10)
    ax.set_ylabel("Number of identified genes", fontsize=10)
    ax.set_title(f"Class 2: Low autocorrelation\n{sig_name}", fontsize=10)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(0, n_genes + 1)

    plt.tight_layout()

    fname = out_dir / f"sig_autocor_rankProd_{sig_name}.pdf"
    fig.savefig(fname, format="pdf", bbox_inches="tight")
    plt.close(fig)

    return fname


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
                except (ValueError, np.linalg.LinAlgError) as e:
                    # Individual heatmap can fail on degenerate data; skip it
                    # but keep building the rest of the PDF. Logged, not silent.
                    logger.warning(
                        "clustermap failed for sig=%r ds=%r: %s: %s",
                        sig, ds, type(e).__name__, e,
                    )
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

    # --- 3. RankProd plots (only if >1 dataset and tables exist) ---
    rank_product_tables = compact_result.get("rank_product_tables", {})
    if n_ds > 1 and rank_product_tables:
        for sig in names_sigs:
            rp_table = rank_product_tables.get(sig)
            if rp_table is not None:
                rp_path = _plot_rank_product(rp_table, sig, out_dir)
                if rp_path:
                    paths.append(rp_path)

    return paths
