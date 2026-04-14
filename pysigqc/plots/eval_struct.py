"""Structural clustering and biclustering heatmaps.

Produces per (dataset, signature):
- Expression clustering heatmaps with hierarchical clustering on genes
- Biclustering heatmaps (continuous z-scores) when biclusters are found
- Binarized biclustering heatmaps
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from ._colors import HEATMAP_CMAP
from ._style import save_pdf, dynamic_fontsize


def plot_struct(
    struct_result: dict,
    names_sigs: list[str],
    names_datasets: list[str],
    out_dir: str | Path,
) -> list[Path]:
    """Generate structural heatmap plots. Returns list of PDF paths."""
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
            try:
                g = sns.clustermap(
                    mat.astype(float),
                    cmap=HEATMAP_CMAP,
                    row_cluster=True, col_cluster=False,
                    figsize=(max(6, mat.shape[1] * 0.15), max(4, mat.shape[0] * 0.3)),
                    xticklabels=False,
                    yticklabels=True if mat.shape[0] <= 50 else False,
                    dendrogram_ratio=(0.15, 0),
                )
                g.ax_heatmap.set_title(f"Expression clustering\n{ds} {sig}", fontsize=9)
                fname = f"sig_eval_struct_clustering_{ds}_{sig}.pdf"
                g.savefig(out_dir / fname, format="pdf", bbox_inches="tight")
                plt.close(g.fig)
                paths.append(out_dir / fname)
            except Exception:
                pass

    # --- 2 & 3. Biclustering heatmaps (continuous + binarized) ---
    any_biclusters = struct_result.get("any_biclusters", False)
    if not any_biclusters:
        return paths

    continuous_figs = []
    binarized_figs = []
    for sig in names_sigs:
        for ds in names_datasets:
            bc = biclust.get(sig, {}).get(ds)
            if bc is None:
                continue
            n_bc = bc.get("n_biclusters", 0)
            if n_bc == 0:
                continue

            # Continuous z-score heatmap
            z = bc.get("z_scores")
            if z is not None and z.shape[0] > 1 and z.shape[1] > 1:
                fig_c, ax_c = plt.subplots(figsize=(max(5, z.shape[1] * 0.1), max(3, z.shape[0] * 0.25)))
                sns.heatmap(z.astype(float), ax=ax_c, cmap=HEATMAP_CMAP,
                            xticklabels=False,
                            yticklabels=True if z.shape[0] <= 50 else False)
                ax_c.set_title(f"Biclustering (z-scores)\n{ds} {sig}", fontsize=9)
                continuous_figs.append(fig_c)

            # Binarized heatmap
            binarized = bc.get("binarized")
            if binarized is not None and binarized.shape[0] > 1 and binarized.shape[1] > 1:
                fig_b, ax_b = plt.subplots(figsize=(max(5, binarized.shape[1] * 0.1), max(3, binarized.shape[0] * 0.25)))
                sns.heatmap(binarized.astype(float), ax=ax_b, cmap="Greys",
                            vmin=0, vmax=1,
                            xticklabels=False,
                            yticklabels=False)
                ax_b.set_title(f"Binarized biclustering\n{ds} {sig}", fontsize=9)
                binarized_figs.append(fig_b)

    if continuous_figs:
        from matplotlib.backends.backend_pdf import PdfPages
        pdf_path = out_dir / "sig_eval_bivariate_clustering.pdf"
        with PdfPages(pdf_path) as pp:
            for f in continuous_figs:
                f.tight_layout()
                pp.savefig(f, bbox_inches="tight")
                plt.close(f)
        paths.append(pdf_path)

    if binarized_figs:
        from matplotlib.backends.backend_pdf import PdfPages
        pdf_path = out_dir / "sig_eval_bivariate_clustering_binarized_maps.pdf"
        with PdfPages(pdf_path) as pp:
            for f in binarized_figs:
                f.tight_layout()
                pp.savefig(f, bbox_inches="tight")
                plt.close(f)
        paths.append(pdf_path)

    return paths
