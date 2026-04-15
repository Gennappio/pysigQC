"""Negative and permutation control boxplots matching R's boxplot.matrix2.R.

R style:
- boxplot with outline=T, outpch=NA (suppress outlier dots)
- stripchart overlay: method='jitter', pch=2 (triangle), cex=0.5
- Dataset colors from dataset_colors()
- Legend: bottom, horizontal, bty='n', cex=0.6
- Adaptive PDF sizing: <=25 metrics -> 7", <=50 -> 10", etc.
- X-axis labels perpendicular (R las=2)
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ._colors import SIG_GENE_COLOR
from ._style import save_pdf


def _boxplot_with_overlay(
    metrics_table,
    original_values: dict,
    title: str,
    out_path: Path,
) -> Path:
    cols = list(metrics_table.columns)
    n = len(cols)

    # R's adaptive PDF sizing
    if n <= 25:
        width = 7
    elif n <= 50:
        width = 10
    elif n <= 75:
        width = 14
    else:
        width = 18

    fig, ax = plt.subplots(figsize=(width, 7))

    # Boxplot (R: outline=T, outpch=NA)
    bp_data = [metrics_table[c].dropna().values for c in cols]
    bp = ax.boxplot(
        bp_data, positions=range(1, n + 1), widths=0.6,
        patch_artist=True, showfliers=False,  # outpch=NA
        medianprops=dict(color="black", linewidth=1),
        whiskerprops=dict(color="black", linewidth=0.8),
        capprops=dict(color="black", linewidth=0.8),
        boxprops=dict(linewidth=0.8),
    )
    for patch in bp["boxes"]:
        patch.set_facecolor("white")
        patch.set_edgecolor("black")

    # Stripchart overlay: jittered triangles (R: pch=2, cex=0.5)
    rng = np.random.default_rng(42)
    obs_vals = [original_values.get(c, np.nan) for c in cols]
    for i, c in enumerate(cols):
        vals = metrics_table[c].dropna().values
        if len(vals) > 0:
            jitter = rng.uniform(-0.15, 0.15, size=len(vals))
            ax.scatter(np.full(len(vals), i + 1) + jitter, vals,
                       marker="^", s=12, c="#808080", edgecolors="none",
                       alpha=0.5, zorder=2)

    # Original signature values: red diamonds (prominent overlay)
    ax.scatter(range(1, n + 1), obs_vals,
               c=SIG_GENE_COLOR, marker="D", s=50, zorder=4,
               edgecolors="black", linewidths=0.5,
               label="Original signature")

    # X-axis labels perpendicular (R: las=2)
    ax.set_xticks(range(1, n + 1))
    max_name = max(len(c) for c in cols) if cols else 10
    label_fs = max(min(6, 7 * 12 / max_name), 4)
    ax.set_xticklabels(cols, rotation=90, fontsize=label_fs)
    ax.set_ylim(-0.05, 1.05)
    ax.set_ylabel("Metric value", fontsize=9)
    ax.set_title(title, fontsize=10)

    # Legend at bottom, horizontal (R: x="bottom", horiz=TRUE, bty="n")
    ax.legend(fontsize=7, loc="lower center", bbox_to_anchor=(0.5, -0.15),
              ncol=2, frameon=False, markerscale=0.8)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="pdf", bbox_inches="tight")
    plt.close(fig)
    return out_path


def plot_negative_control(
    neg_result: dict,
    original_radar_result: dict,
    names_sigs: list[str],
    names_datasets: list[str],
    out_dir: str | Path,
) -> list[Path]:
    out_dir = Path(out_dir)
    paths: list[Path] = []

    output_table = original_radar_result.get("output_table")

    for control_type in ("negative_controls", "permutation_controls"):
        controls = neg_result.get(control_type, {})
        subdir = "negative_control" if "negative" in control_type else "permutation_control"
        for ds in names_datasets:
            for sig in names_sigs:
                data = controls.get(ds, {}).get(sig)
                if data is None:
                    continue
                metrics_table = data.get("metrics_table")
                if metrics_table is None or metrics_table.empty:
                    continue

                orig = {}
                if output_table is not None:
                    row_key = f"{sig}_{ds}"
                    if row_key in output_table.index:
                        orig = output_table.loc[row_key].to_dict()

                pdf_path = out_dir / subdir / ds / sig / "boxplot_metrics.pdf"
                paths.append(
                    _boxplot_with_overlay(
                        metrics_table, orig,
                        title=f"{subdir.replace('_', ' ').title()}\n{ds} \u2014 {sig}",
                        out_path=pdf_path,
                    )
                )

    return paths
