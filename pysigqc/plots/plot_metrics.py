"""Score comparison plots matching R's compare_metrics_loc.R styling.

Per signature produces ``sig_compare_metrics_<sig>.pdf`` with layout
mfcol=c(4, n_datasets):
  Row 1: Median vs Mean (smoothScatter + red density)
  Row 2: Mean vs PCA1
  Row 3: PCA1 vs Median
  Row 4: PCA scree barplot

Also produces:
- ``sig_qq_plots_<sig>.pdf``: QQ plots with layout mfcol=c(3, n_datasets)
- ``sig_compare_ES_metrics_<sig>.pdf``: Enrichment score comparisons
  (GSVA vs ssGSEA, ssGSEA vs PLAGE, PLAGE vs GSVA)
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import gaussian_kde, spearmanr, probplot

from ._style import save_pdf, dynamic_fontsize

# R's jet colormap
_JET_COLORS = ["#00007F", "#0000FF", "#007FFF", "#00FFFF",
               "#7FFF7F", "#FFFF00", "#FF7F00", "#FF0000", "#7F0000"]
_JET_CMAP = LinearSegmentedColormap.from_list("r_jet", _JET_COLORS, N=256)


def _smoothscatter(ax, x, y, xlabel, ylabel, title, font):
    """R-style smoothScatter: 2D density heatmap + black points + red density line.

    Matches R's graphics::smoothScatter() which:
    1. Creates a filled 2D kernel density image (blue-cyan-yellow-red)
    2. Overlays individual points as small black dots
    3. Overlays 1D density of x-axis as red line on secondary y-axis
    """
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]

    if len(x) < 3:
        ax.text(0.5, 0.5, "Insufficient data", transform=ax.transAxes,
                ha="center", va="center", fontsize=9, color="grey")
        ax.set_title(title, fontsize=max(font, 10))
        return

    # Compute axis limits with small padding
    x_min, x_max = x.min(), x.max()
    y_min, y_max = y.min(), y.max()
    x_pad = (x_max - x_min) * 0.05 or 0.1
    y_pad = (y_max - y_min) * 0.05 or 0.1

    # Create 2D density heatmap (filled image like R's smoothScatter)
    try:
        # Create a grid for the 2D KDE
        n_grid = 100
        xx = np.linspace(x_min - x_pad, x_max + x_pad, n_grid)
        yy = np.linspace(y_min - y_pad, y_max + y_pad, n_grid)
        XX, YY = np.meshgrid(xx, yy)
        positions = np.vstack([XX.ravel(), YY.ravel()])

        # Compute 2D KDE
        xy = np.vstack([x, y])
        kde = gaussian_kde(xy)
        Z = kde(positions).reshape(XX.shape)

        # Plot as filled contour/image (matching R's blue-to-yellow-to-red)
        ax.imshow(Z, origin='lower', aspect='auto',
                  extent=[x_min - x_pad, x_max + x_pad, y_min - y_pad, y_max + y_pad],
                  cmap=_JET_CMAP, interpolation='bilinear')

        # Overlay individual points as small black dots (R uses pch='.')
        ax.scatter(x, y, c='black', s=3, marker='.', edgecolors='none', alpha=0.8)

    except np.linalg.LinAlgError:
        # Fallback: just scatter plot with density coloring
        ax.scatter(x, y, s=4, c="blue", edgecolors="none", rasterized=True)

    # Red density overlay on right axis (1D density of x)
    # Use smaller bandwidth to match R's density() which shows more detail/ripples
    try:
        kde_1d = gaussian_kde(x)
        # Reduce bandwidth by 50% to show more ripples like R's density()
        kde_1d.set_bandwidth(kde_1d.factor * 0.5)
        xs = np.linspace(x_min, x_max, 300)
        dens = kde_1d(xs)
        ax2 = ax.twinx()
        ax2.plot(xs, dens, color="red", lw=2, alpha=0.9)
        ax2.set_ylabel("Density", fontsize=8, color="red")
        ax2.tick_params(axis="y", labelsize=6, colors="red")
        ax2.set_ylim(0, dens.max() * 1.1)
    except np.linalg.LinAlgError:
        pass

    # Compute and display Spearman correlation
    rho, _ = spearmanr(x, y)

    # Title with rho annotation (R puts rho in top-right corner)
    ax.set_title(title, fontsize=max(font, 10), fontweight='bold')
    ax.annotate(f"rho = {rho:.1f}", xy=(0.95, 0.95), xycoords='axes fraction',
                ha='right', va='top', fontsize=8,
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))

    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.tick_params(labelsize=7)


def plot_metrics(
    metrics_result: dict,
    names_sigs: list[str],
    names_datasets: list[str],
    out_dir: str | Path,
) -> list[Path]:
    scores = metrics_result["scores"]
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    n_ds = len(names_datasets)
    max_label = max(len(f"{ds} {sig}") for sig in names_sigs for ds in names_datasets)
    # Increase base font size - R uses larger titles
    font = max(9.0, min(11.0, 6 * 10 / max_label))

    paths: list[Path] = []

    for sig in names_sigs:
        # --- Main comparison PDF: mfcol=c(4, n_ds) ---
        fig, axes = plt.subplots(4, n_ds, figsize=(3 * n_ds, 10), squeeze=False)

        for di, ds in enumerate(names_datasets):
            sc = scores.get(sig, {}).get(ds, {})
            med = np.asarray(sc.get("med_scores", []), dtype=float)
            mean = np.asarray(sc.get("mean_scores", []), dtype=float)
            pca1 = np.asarray(sc.get("pca1_scores", []), dtype=float) if sc.get("pca1_scores") is not None else np.array([])
            props = sc.get("props_of_variances")

            # Row 0: Median vs Mean
            _smoothscatter(axes[0, di], med, mean, "Median", "Mean",
                           f"Median vs Mean\n{ds}", font)
            # Row 1: Mean vs PCA1
            if len(pca1) > 0:
                _smoothscatter(axes[1, di], mean, pca1, "Mean", "PCA1",
                               f"Mean vs PCA1\n{ds}", font)
            else:
                axes[1, di].text(0.5, 0.5, "No PCA", transform=axes[1, di].transAxes,
                                 ha="center", va="center", fontsize=7, color="grey")
                axes[1, di].set_title(f"Mean vs PCA1\n{ds}", fontsize=font)
            # Row 2: PCA1 vs Median
            if len(pca1) > 0:
                _smoothscatter(axes[2, di], pca1, med, "PCA1", "Median",
                               f"PCA1 vs Median\n{ds}", font)
            else:
                axes[2, di].text(0.5, 0.5, "No PCA", transform=axes[2, di].transAxes,
                                 ha="center", va="center", fontsize=7, color="grey")
                axes[2, di].set_title(f"PCA1 vs Median\n{ds}", fontsize=font)
            # Row 3: Scree plot (R: barplot, default style)
            if props is not None and len(props) > 0:
                n_show = min(10, len(props))
                axes[3, di].bar(range(1, n_show + 1), props[:n_show],
                                color="#D3D3D3", edgecolor="black", linewidth=0.5)
                axes[3, di].set_xlabel("PC", fontsize=7)
                axes[3, di].set_ylabel("Prop. variance", fontsize=7)
                axes[3, di].set_title(f"PCA vs proportion\nof variance\n{ds}", fontsize=font)
            else:
                axes[3, di].text(0.5, 0.5, "No PCA", transform=axes[3, di].transAxes,
                                 ha="center", va="center", fontsize=7, color="grey")
                axes[3, di].set_title(f"Scree plot\n{ds}", fontsize=font)

        fig.suptitle(sig, fontsize=11, y=1.01)
        fig.tight_layout()
        paths.append(save_pdf(fig, out_dir, f"sig_compare_metrics_{sig}.pdf"))

    # --- QQ plots: mfcol=c(3, n_ds) ---
    for sig in names_sigs:
        fig, axes = plt.subplots(3, n_ds, figsize=(3 * n_ds, 9), squeeze=False)
        for di, ds in enumerate(names_datasets):
            sc = scores.get(sig, {}).get(ds, {})
            for ri, (label, key) in enumerate(
                [("Median score", "med_scores"),
                 ("Mean score", "mean_scores"),
                 ("PCA1 score", "pca1_scores")]
            ):
                ax = axes[ri, di]
                vals = sc.get(key)
                if vals is not None:
                    vals = np.asarray(vals, dtype=float)
                    vals = vals[np.isfinite(vals)]
                    if len(vals) > 2:
                        probplot(vals, dist="norm", plot=ax)
                        ax.set_title(f"{label}\n{ds}", fontsize=font)
                        # Style the QQ line
                        for line in ax.get_lines():
                            if line.get_linestyle() == '--':
                                line.set_color("red")
                                line.set_linewidth(1)
                        continue
                ax.text(0.5, 0.5, "N/A", transform=ax.transAxes,
                        ha="center", va="center", fontsize=7, color="grey")
                ax.set_title(f"{label}\n{ds}", fontsize=font)

        fig.suptitle(sig, fontsize=11, y=1.01)
        fig.tight_layout()
        paths.append(save_pdf(fig, out_dir, f"sig_qq_plots_{sig}.pdf"))

    # --- Scoring metrics correlation heatmaps ---
    score_cor_mats = metrics_result.get("score_cor_mats", {})
    if score_cor_mats:
        import seaborn as sns
        for sig in names_sigs:
            for ds in names_datasets:
                key = f"{ds}_{sig}"
                cor_mat = score_cor_mats.get(key)
                if cor_mat is None or cor_mat.shape[0] < 2:
                    continue
                fig, ax = plt.subplots(figsize=(4, 4))
                sns.heatmap(
                    cor_mat, vmin=-1, vmax=1, cmap="RdBu_r", center=0,
                    annot=True, fmt=".2f", square=True, ax=ax,
                    linewidths=0.5, cbar_kws={"shrink": 0.8, "label": "Spearman rho"},
                )
                ax.set_title(f"Scoring Metric Correlation\n{ds}  {sig}", fontsize=9)
                fig.tight_layout()
                paths.append(save_pdf(fig, out_dir, f"scoring_metrics_corr_{ds}_{sig}.pdf"))

    # --- Gaussian Mixture Model BIC plots ---
    mixture_models = metrics_result.get("mixture_models", {})
    if mixture_models:
        for sig in names_sigs:
            fig, axes = plt.subplots(3, n_ds, figsize=(3 * n_ds, 9), squeeze=False)
            has_data = False
            for di, ds in enumerate(names_datasets):
                mm = mixture_models.get(sig, {}).get(ds, {})
                for ri, score_name in enumerate(["median", "mean", "pca1"]):
                    ax = axes[ri, di]
                    entry = mm.get(score_name)
                    bic_values = None
                    if isinstance(entry, dict):
                        bic_values = entry.get("bic_values", [])
                    if bic_values:
                        ks, bics = zip(*bic_values)
                        ax.plot(ks, bics, "o-", color="black", markersize=4, linewidth=1)
                        ax.set_xlabel("Number of components", fontsize=7)
                        ax.set_ylabel("BIC", fontsize=7)
                        ax.set_title(f"{score_name.title()} score\n{ds}", fontsize=font)
                        ax.tick_params(labelsize=6)
                        has_data = True
                    else:
                        ax.text(0.5, 0.5, "N/A", transform=ax.transAxes,
                                ha="center", va="center", fontsize=7, color="grey")
                        ax.set_title(f"{score_name.title()} score\n{ds}", fontsize=font)

            if has_data:
                fig.suptitle(f"Gaussian Mixture Model BIC\n{sig}", fontsize=11, y=1.02)
                fig.tight_layout()
                paths.append(save_pdf(fig, out_dir, f"sig_gaussian_mixture_model_{sig}.pdf"))
            else:
                plt.close(fig)

    # --- Enrichment Score Comparison plots: sig_compare_ES_metrics_<sig>.pdf ---
    # Layout: 3 rows (GSVA vs ssGSEA, ssGSEA vs PLAGE, PLAGE vs GSVA) x n_datasets
    enrichment_scores = metrics_result.get("enrichment_scores", {})
    if enrichment_scores:
        for sig in names_sigs:
            fig, axes = plt.subplots(3, n_ds, figsize=(3 * n_ds, 7.5), squeeze=False)
            has_data = False

            for di, ds in enumerate(names_datasets):
                es = enrichment_scores.get(sig, {}).get(ds, {})
                gsva = es.get("gsva")
                ssgsea = es.get("ssgsea")
                plage = es.get("plage")

                # Row 0: GSVA vs ssGSEA
                ax = axes[0, di]
                if gsva is not None and ssgsea is not None:
                    gsva_arr = np.asarray(gsva, dtype=float)
                    ssgsea_arr = np.asarray(ssgsea, dtype=float)
                    mask = np.isfinite(gsva_arr) & np.isfinite(ssgsea_arr)
                    if mask.sum() > 2:
                        _smoothscatter(ax, gsva_arr[mask], ssgsea_arr[mask],
                                      "GSVA", "ssGSEA", f"GSVA vs ssGSEA\n{ds}", font)
                        has_data = True
                    else:
                        ax.text(0.5, 0.5, "Insufficient data", transform=ax.transAxes,
                                ha="center", va="center", fontsize=7, color="grey")
                        ax.set_title(f"GSVA vs ssGSEA\n{ds}", fontsize=font)
                else:
                    ax.text(0.5, 0.5, "N/A", transform=ax.transAxes,
                            ha="center", va="center", fontsize=7, color="grey")
                    ax.set_title(f"GSVA vs ssGSEA\n{ds}", fontsize=font)

                # Row 1: ssGSEA vs PLAGE
                ax = axes[1, di]
                if ssgsea is not None and plage is not None:
                    ssgsea_arr = np.asarray(ssgsea, dtype=float)
                    plage_arr = np.asarray(plage, dtype=float)
                    mask = np.isfinite(ssgsea_arr) & np.isfinite(plage_arr)
                    if mask.sum() > 2:
                        _smoothscatter(ax, ssgsea_arr[mask], plage_arr[mask],
                                      "ssGSEA", "PLAGE", f"ssGSEA vs PLAGE\n{ds}", font)
                        has_data = True
                    else:
                        ax.text(0.5, 0.5, "Insufficient data", transform=ax.transAxes,
                                ha="center", va="center", fontsize=7, color="grey")
                        ax.set_title(f"ssGSEA vs PLAGE\n{ds}", fontsize=font)
                else:
                    ax.text(0.5, 0.5, "N/A", transform=ax.transAxes,
                            ha="center", va="center", fontsize=7, color="grey")
                    ax.set_title(f"ssGSEA vs PLAGE\n{ds}", fontsize=font)

                # Row 2: PLAGE vs GSVA
                ax = axes[2, di]
                if plage is not None and gsva is not None:
                    plage_arr = np.asarray(plage, dtype=float)
                    gsva_arr = np.asarray(gsva, dtype=float)
                    mask = np.isfinite(plage_arr) & np.isfinite(gsva_arr)
                    if mask.sum() > 2:
                        _smoothscatter(ax, plage_arr[mask], gsva_arr[mask],
                                      "PLAGE", "GSVA", f"PLAGE vs GSVA\n{ds}", font)
                        has_data = True
                    else:
                        ax.text(0.5, 0.5, "Insufficient data", transform=ax.transAxes,
                                ha="center", va="center", fontsize=7, color="grey")
                        ax.set_title(f"PLAGE vs GSVA\n{ds}", fontsize=font)
                else:
                    ax.text(0.5, 0.5, "N/A", transform=ax.transAxes,
                            ha="center", va="center", fontsize=7, color="grey")
                    ax.set_title(f"PLAGE vs GSVA\n{ds}", fontsize=font)

            if has_data:
                fig.suptitle(f"Enrichment Score Comparison\n{sig}", fontsize=11, y=1.02)
                fig.tight_layout()
                paths.append(save_pdf(fig, out_dir, f"sig_compare_ES_metrics_{sig}.pdf"))
            else:
                plt.close(fig)

    return paths
