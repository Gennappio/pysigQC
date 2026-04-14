"""Score comparison plots: scatter pairs (mean/median/PCA1), scree plot,
Gaussian mixture BIC plots, QQ plots.

Produces one multi-page PDF per signature: ``sig_compare_metrics_<sig>.pdf``
containing 4 rows of subplots (one column per dataset):
  Row 1: Median vs Mean scatter
  Row 2: Mean vs PCA1 scatter
  Row 3: PCA1 vs Median scatter
  Row 4: PCA scree plot

Separate PDFs for GMM BIC and QQ plots.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import spearmanr, probplot

from ._colors import JET_CMAP
from ._style import save_pdf, dynamic_fontsize


def _scatter_with_rho(ax, x, y, xlabel, ylabel, title, font):
    """Density-colored scatter with Spearman rho annotation."""
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) > 2:
        ax.hexbin(x, y, gridsize=30, cmap=JET_CMAP, mincnt=1)
        rho, _ = spearmanr(x, y)
        ax.text(0.05, 0.95, f"rho={rho:.3f}", transform=ax.transAxes,
                fontsize=7, va="top",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7))
    else:
        ax.text(0.5, 0.5, "Insufficient data", transform=ax.transAxes,
                ha="center", va="center", fontsize=8)
    ax.set_xlabel(xlabel, fontsize=7)
    ax.set_ylabel(ylabel, fontsize=7)
    ax.set_title(title, fontsize=font)


def plot_metrics(
    metrics_result: dict,
    names_sigs: list[str],
    names_datasets: list[str],
    out_dir: str | Path,
) -> list[Path]:
    """Generate score comparison plots. Returns list of PDF paths."""
    scores = metrics_result["scores"]
    pca_results = metrics_result.get("pca_results", {})
    mixture = metrics_result.get("mixture_models", {})
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    n_ds = len(names_datasets)
    max_label = max(len(f"{ds} {sig}") for sig in names_sigs for ds in names_datasets)
    font = dynamic_fontsize(max_label)

    paths: list[Path] = []

    for sig in names_sigs:
        # --- Main comparison PDF: 4 rows x n_ds columns ---
        fig, axes = plt.subplots(4, n_ds, figsize=(3 * n_ds, 10), squeeze=False)

        for di, ds in enumerate(names_datasets):
            sc = scores.get(sig, {}).get(ds, {})
            med = np.asarray(sc.get("med_scores", []), dtype=float)
            mean = np.asarray(sc.get("mean_scores", []), dtype=float)
            pca1 = np.asarray(sc.get("pca1_scores", []), dtype=float) if sc.get("pca1_scores") is not None else np.array([])
            props = sc.get("props_of_variances")

            # Row 0: Median vs Mean
            _scatter_with_rho(axes[0, di], med, mean, "Median", "Mean",
                              f"Med vs Mean\n{ds}", font)
            # Row 1: Mean vs PCA1
            if len(pca1) > 0:
                _scatter_with_rho(axes[1, di], mean, pca1, "Mean", "PCA1",
                                  f"Mean vs PCA1\n{ds}", font)
            else:
                axes[1, di].text(0.5, 0.5, "No PCA", transform=axes[1, di].transAxes,
                                 ha="center", va="center")
                axes[1, di].set_title(f"Mean vs PCA1\n{ds}", fontsize=font)
            # Row 2: PCA1 vs Median
            if len(pca1) > 0:
                _scatter_with_rho(axes[2, di], pca1, med, "PCA1", "Median",
                                  f"PCA1 vs Median\n{ds}", font)
            else:
                axes[2, di].text(0.5, 0.5, "No PCA", transform=axes[2, di].transAxes,
                                 ha="center", va="center")
                axes[2, di].set_title(f"PCA1 vs Median\n{ds}", fontsize=font)
            # Row 3: Scree plot
            if props is not None and len(props) > 0:
                n_show = min(10, len(props))
                axes[3, di].bar(range(1, n_show + 1), props[:n_show], color="#4682B4")
                axes[3, di].set_xlabel("PC", fontsize=7)
                axes[3, di].set_ylabel("Prop. variance", fontsize=7)
                axes[3, di].set_title(f"Scree plot\n{ds}", fontsize=font)
            else:
                axes[3, di].text(0.5, 0.5, "No PCA", transform=axes[3, di].transAxes,
                                 ha="center", va="center")
                axes[3, di].set_title(f"Scree plot\n{ds}", fontsize=font)

        fig.suptitle(sig, fontsize=12, y=1.01)
        fig.tight_layout()
        paths.append(save_pdf(fig, out_dir, f"sig_compare_metrics_{sig}.pdf"))

    # --- GMM BIC plots ---
    if mixture:
        for sig in names_sigs:
            figs_bic = []
            for ds in names_datasets:
                mm = mixture.get(sig, {}).get(ds, {})
                fig_b, bic_axes = plt.subplots(1, 3, figsize=(9, 3))
                for idx, (label, model) in enumerate(
                    [("Median", mm.get("median")),
                     ("Mean", mm.get("mean")),
                     ("PCA1", mm.get("pca1"))]
                ):
                    ax = bic_axes[idx]
                    if model is not None and hasattr(model, "bic_"):
                        # sklearn GaussianMixture doesn't store per-k BIC;
                        # just show the selected model's n_components
                        ax.text(0.5, 0.5,
                                f"K={model.n_components}",
                                transform=ax.transAxes, ha="center", va="center",
                                fontsize=14)
                    else:
                        ax.text(0.5, 0.5, "N/A", transform=ax.transAxes,
                                ha="center", va="center")
                    ax.set_title(f"{label}\n{ds}", fontsize=8)
                fig_b.suptitle(f"GMM — {sig}", fontsize=10)
                fig_b.tight_layout()
                figs_bic.append(fig_b)
            if figs_bic:
                pdf_path = out_dir / f"sig_gaussian_mixture_model_{sig}.pdf"
                with PdfPages(pdf_path) as pp:
                    for f in figs_bic:
                        pp.savefig(f, bbox_inches="tight")
                        plt.close(f)
                paths.append(pdf_path)

    # --- QQ plots ---
    for sig in names_sigs:
        fig_qq, qq_axes = plt.subplots(3, n_ds, figsize=(3 * n_ds, 9), squeeze=False)
        for di, ds in enumerate(names_datasets):
            sc = scores.get(sig, {}).get(ds, {})
            for ri, (label, key) in enumerate(
                [("Median", "med_scores"), ("Mean", "mean_scores"), ("PCA1", "pca1_scores")]
            ):
                ax = qq_axes[ri, di]
                vals = sc.get(key)
                if vals is not None:
                    vals = np.asarray(vals, dtype=float)
                    vals = vals[np.isfinite(vals)]
                    if len(vals) > 2:
                        probplot(vals, dist="norm", plot=ax)
                        ax.set_title(f"QQ {label}\n{ds}", fontsize=font)
                        continue
                ax.text(0.5, 0.5, "N/A", transform=ax.transAxes,
                        ha="center", va="center")
                ax.set_title(f"QQ {label}\n{ds}", fontsize=font)
        fig_qq.suptitle(sig, fontsize=12, y=1.01)
        fig_qq.tight_layout()
        paths.append(save_pdf(fig_qq, out_dir, f"sig_qq_plots_{sig}.pdf"))

    return paths
