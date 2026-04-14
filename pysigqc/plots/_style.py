"""Shared plotting helpers: figure layout, PDF saving, font scaling."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # non-interactive backend for PDF generation
import matplotlib.pyplot as plt
import numpy as np


def figure_grid(
    n_rows: int,
    n_cols: int,
    cell_w: float = 4.0,
    cell_h: float = 4.0,
) -> tuple[plt.Figure, np.ndarray]:
    """Create a figure with an *n_rows* x *n_cols* subplot grid.

    Returns ``(fig, axes)`` where *axes* is always 2-D.
    """
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(cell_w * n_cols, cell_h * n_rows),
        squeeze=False,
    )
    return fig, axes


def save_pdf(fig: plt.Figure, out_dir: str | Path, filename: str) -> Path:
    """Save *fig* as a PDF in *out_dir*, close it, and return the path."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / filename
    fig.savefig(path, format="pdf", bbox_inches="tight")
    plt.close(fig)
    return path


def dynamic_fontsize(max_label_length: int, base: float = 10.0) -> float:
    """Compute a font size that shrinks for long labels, matching R's heuristic."""
    if max_label_length <= 0:
        return base
    return min(base, 4 * 10 / max_label_length)
