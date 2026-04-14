"""Radar (spider) chart summarising all 14 QC metrics.

Produces ``sig_radarplot.pdf`` — one overlaid polygon per (sig, dataset)
combination, colored by dataset and styled by signature.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ._colors import dataset_colors, sig_linestyle
from ._style import save_pdf


# Metric labels for the 14 spokes
_METRIC_LABELS = [
    "Rel. med. SD", "Abs. skewness", "Top 10% CV", "Top 25% CV",
    "Top 50% CV", "CV ratio", "Med. prop. non-NA", "Med. prop. expr.",
    "Autocorrelation", "rho(Mean,Med)", "rho(PCA1,Med)", "rho(Mean,PCA1)",
    "Prop. var. PC1", "Std. robustness",
]


def plot_radar(
    radar_result: dict,
    names_sigs: list[str],
    names_datasets: list[str],
    out_dir: str | Path,
) -> Path:
    """Generate summary radar chart. Returns PDF path."""
    # radar_plot_mat has rows: [max, min, ...data rows]
    radar_plot_mat = np.asarray(radar_result["radar_plot_mat"], dtype=float)
    areas = radar_result.get("areas", [])
    legend_labels = radar_result.get("legend_labels", [])

    n_ds = len(names_datasets)
    n_sigs = len(names_sigs)
    n_metrics = radar_plot_mat.shape[1]
    ds_colors = dataset_colors(n_ds)

    # Angular positions (start at top, go clockwise)
    angles = np.linspace(0, 2 * np.pi, n_metrics, endpoint=False).tolist()
    angles += angles[:1]  # close the polygon

    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(polar=True))
    ax.set_theta_offset(np.pi / 2)   # start at 12 o'clock
    ax.set_theta_direction(-1)        # clockwise

    # Draw grid at 0.2 intervals
    ax.set_ylim(0, 1)
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(["0.2", "0.4", "0.6", "0.8", "1.0"], fontsize=6, color="grey")

    # Spoke labels
    labels = _METRIC_LABELS[:n_metrics]
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels, fontsize=7)

    # Plot each (sig, dataset) polygon
    data_rows = radar_plot_mat[2:]  # skip max/min rows
    row_idx = 0
    handles = []
    labs = []
    for ki, sig in enumerate(names_sigs):
        for di, ds in enumerate(names_datasets):
            if row_idx >= len(data_rows):
                break
            vals = np.abs(data_rows[row_idx])  # abs() matches R behavior
            vals_closed = np.concatenate([vals, [vals[0]]])  # close polygon

            ls = sig_linestyle(ki)
            line, = ax.plot(angles, vals_closed, color=ds_colors[di],
                            linestyle=ls, lw=1.8, alpha=0.85)
            ax.fill(angles, vals_closed, color=ds_colors[di], alpha=0.06)

            label = legend_labels[row_idx] if row_idx < len(legend_labels) else f"{ds} {sig}"
            handles.append(line)
            labs.append(label)
            row_idx += 1

    # Legend outside the chart
    ax.legend(handles, labs, loc="upper left", bbox_to_anchor=(1.15, 1.05),
              fontsize=7, frameon=False)

    fig.tight_layout()
    return save_pdf(fig, out_dir, "sig_radarplot.pdf")
