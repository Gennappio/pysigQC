"""Negative and permutation control boxplots.

Produces per (dataset, signature):
- ``negative_control/<ds>/<sig>/boxplot_metrics.pdf``
- ``permutation_control/<ds>/<sig>/boxplot_metrics.pdf``

Each shows boxplots of the null distribution (from resampling) with the
original signature's metric overlaid as a red diamond.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ._colors import SIG_GENE_COLOR, dataset_colors
from ._style import save_pdf


def _boxplot_with_overlay(
    metrics_table,          # DataFrame: n_resamples x n_metrics
    original_values: dict,  # metric_name -> float (observed value)
    title: str,
    out_path: Path,
) -> Path:
    """Draw boxplots of null distribution with observed values overlaid."""
    cols = list(metrics_table.columns)
    n = len(cols)
    width = max(7, n * 0.5)
    fig, ax = plt.subplots(figsize=(width, 5))

    # Boxplot
    bp_data = [metrics_table[c].dropna().values for c in cols]
    bp = ax.boxplot(bp_data, positions=range(n), widths=0.6,
                    patch_artist=True, zorder=1)
    for patch in bp["boxes"]:
        patch.set_facecolor("#D3D3D3")
        patch.set_edgecolor("#808080")

    # Overlay observed values
    obs = [original_values.get(c, np.nan) for c in cols]
    ax.scatter(range(n), obs, c=SIG_GENE_COLOR, marker="D", s=40, zorder=3,
               edgecolors="black", linewidths=0.5, label="Original signature")

    ax.set_xticks(range(n))
    ax.set_xticklabels(cols, rotation=90, fontsize=6)
    ax.set_ylim(-0.05, 1.05)
    ax.set_ylabel("Metric value")
    ax.set_title(title, fontsize=9)
    ax.legend(fontsize=7, loc="upper right")

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
    """Generate negative/permutation control boxplots. Returns PDF paths."""
    out_dir = Path(out_dir)
    paths: list[Path] = []

    # Get original radar values from the output_table
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

                # Extract original values for overlay
                orig = {}
                if output_table is not None:
                    row_key = f"{sig}_{ds}"
                    if row_key in output_table.index:
                        orig = output_table.loc[row_key].to_dict()

                pdf_path = out_dir / subdir / ds / sig / "boxplot_metrics.pdf"
                paths.append(
                    _boxplot_with_overlay(
                        metrics_table, orig,
                        title=f"{subdir.replace('_', ' ').title()}\n{ds} — {sig}",
                        out_path=pdf_path,
                    )
                )

    return paths
