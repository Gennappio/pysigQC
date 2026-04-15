"""Structural clustering and biclustering heatmaps matching R's eval_struct_loc.R.

Expression heatmaps: seaborn clustermap with row dendrograms, no column
clustering (matching ComplexHeatmap with show_column_dend=F).

Biclustering: blue-white-red continuous heatmap, and binary black/white
heatmap (matching gplots::heatmap.2 and biclust::heatmapBC).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

from ._style import save_pdf


_HEATMAP_CMAP = "RdBu_r"


def plot_struct(
    struct_result: dict,
    names_sigs: list[str],
    names_datasets: list[str],
    out_dir: str | Path,
) -> list[Path]:
    sig_scores = struct_result.get("sig_scores_all_mats", {})
    biclust = struct_result.get("biclust_results", {})
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    paths: list[Path] = []

    # --- 1. Expression clustering heatmaps ---
    for sig in names_sigs:
        for ds in names_datasets:
            mat = sig_scores.get(sig, {}).get(ds)
            if mat is None or mat.shape[0] < 2 or mat.shape[1] < 2:
                continue
            n_genes = mat.shape[0]
            n_samples = mat.shape[1]
            label_fs = max(min(5, 4 * 4 / n_genes * 10), 3)
            fig_w = max(5, n_samples * 0.08 + 2)
            fig_h = max(4, n_genes * 0.25 + 2)

            try:
                g = sns.clustermap(
                    mat.astype(float),
                    cmap=_HEATMAP_CMAP,
                    row_cluster=True, col_cluster=False,
                    figsize=(fig_w, fig_h),
                    xticklabels=False,
                    yticklabels=True if n_genes <= 50 else False,
                    dendrogram_ratio=(0.15, 0),
                    linewidths=0, linecolor="white",
                    cbar_kws={"shrink": 0.5},
                )
                g.ax_heatmap.tick_params(labelsize=label_fs)
                g.ax_heatmap.set_title(f"{ds}  {sig}", fontsize=9, pad=8)
                fname = f"sig_eval_struct_clustering_{ds}_{sig}.pdf"
                g.savefig(out_dir / fname, format="pdf", bbox_inches="tight")
                plt.close(g.fig)
                paths.append(out_dir / fname)
            except Exception:
                pass

    # --- 2 & 3. Biclustering heatmaps ---
    any_biclusters = struct_result.get("any_biclusters", False)
    if not any_biclusters:
        return paths

    continuous_figs = []
    binarized_figs = []
    for sig in names_sigs:
        for ds in names_datasets:
            bc = biclust.get(sig, {}).get(ds)
            if bc is None or bc.get("n_biclusters", 0) == 0:
                continue

            z = bc.get("z_scores")
            if z is not None and z.shape[0] > 1 and z.shape[1] > 1:
                fig_c, ax_c = plt.subplots(
                    figsize=(max(5, z.shape[1] * 0.08 + 1),
                             max(3, z.shape[0] * 0.2 + 1)))
                sns.heatmap(z.astype(float), ax=ax_c, cmap=_HEATMAP_CMAP,
                            center=0, xticklabels=False,
                            yticklabels=True if z.shape[0] <= 50 else False,
                            linewidths=0, cbar_kws={"shrink": 0.5, "label": "Z-score"})
                ax_c.set_title(f"Biclustering (z-scores)\n{ds}  {sig}", fontsize=9)
                continuous_figs.append(fig_c)

            binarized = bc.get("binarized")
            if binarized is not None and binarized.shape[0] > 1 and binarized.shape[1] > 1:
                fig_b, ax_b = plt.subplots(
                    figsize=(max(5, binarized.shape[1] * 0.08 + 1),
                             max(3, binarized.shape[0] * 0.2 + 1)))
                sns.heatmap(binarized.astype(float), ax=ax_b, cmap="Greys",
                            vmin=0, vmax=1, xticklabels=False, yticklabels=False,
                            linewidths=0, cbar=False)
                ax_b.set_title(f"Binarized biclustering\n{ds}  {sig}", fontsize=9)
                binarized_figs.append(fig_b)

    if continuous_figs:
        pdf_path = out_dir / "sig_eval_bivariate_clustering.pdf"
        with PdfPages(pdf_path) as pp:
            for f in continuous_figs:
                f.tight_layout()
                pp.savefig(f, bbox_inches="tight")
                plt.close(f)
        paths.append(pdf_path)

    if binarized_figs:
        pdf_path = out_dir / "sig_eval_bivariate_clustering_binarized_maps.pdf"
        with PdfPages(pdf_path) as pp:
            for f in binarized_figs:
                f.tight_layout()
                pp.savefig(f, bbox_inches="tight")
                plt.close(f)
        paths.append(pdf_path)

    return paths
