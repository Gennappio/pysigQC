"""Decoupled plotting layer for pysigqc.

Usage::

    from pysigqc_joblib.pipeline import run_pipeline
    from pysigqc.plots import plot_all

    result = run_pipeline(sigs, names_sigs, datasets, names_datasets)
    paths = plot_all(result, out_dir="plots/")

Or standalone on previously saved results::

    from pysigqc.plots import plot_all
    paths = plot_all(loaded_dict, out_dir="plots/")
"""

from __future__ import annotations

import logging
from pathlib import Path

from .eval_var import plot_var
from .eval_expr import plot_expr
from .eval_stan import plot_stan
from .eval_compactness import plot_compactness
from .eval_struct import plot_struct
from .compare_metrics import plot_metrics
from .radar import plot_radar
from .negative_control import plot_negative_control

logger = logging.getLogger(__name__)


def plot_all(
    pipeline_result: dict,
    out_dir: str | Path,
    names_sigs: list[str] | None = None,
    names_datasets: list[str] | None = None,
) -> dict[str, list[Path]]:
    """Generate all plots from a pipeline result dict.

    Parameters
    ----------
    pipeline_result : dict
        As returned by ``run_pipeline()`` or assembled from saved files.
    out_dir : str or Path
        Directory for PDF output (created if needed).
    names_sigs, names_datasets : list[str], optional
        If *None*, inferred from ``pipeline_result["radar_values"]`` keys.

    Returns
    -------
    dict mapping module name -> list of generated PDF paths.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Infer names from radar_values if not provided
    if names_sigs is None or names_datasets is None:
        rv = pipeline_result.get("radar_values", {})
        if names_sigs is None:
            names_sigs = list(rv.keys())
        if names_datasets is None and names_sigs:
            names_datasets = list(rv[names_sigs[0]].keys())

    generated: dict[str, list[Path]] = {}

    modules = [
        ("var", plot_var, "var_result"),
        ("expr", plot_expr, "expr_result"),
        ("stan", plot_stan, "stan_result"),
        ("compactness", plot_compactness, "compact_result"),
        ("struct", plot_struct, "struct_result"),
        ("metrics", plot_metrics, "metrics_result"),
        ("radar", plot_radar, "radar_result"),
    ]

    for name, func, key in modules:
        if key not in pipeline_result:
            logger.info("Skipping %s plots — %r not in pipeline_result", name, key)
            continue
        try:
            result = func(
                pipeline_result[key],
                names_sigs=names_sigs,
                names_datasets=names_datasets,
                out_dir=out_dir,
            )
            if isinstance(result, Path):
                result = [result]
            generated[name] = list(result)
        except Exception:
            logger.exception("Error generating %s plots", name)

    # Negative control is optional and uses a different key structure
    neg_key = "negative_result"
    if neg_key in pipeline_result:
        try:
            result = plot_negative_control(
                pipeline_result[neg_key],
                original_radar_result=pipeline_result.get("radar_result", {}),
                names_sigs=names_sigs,
                names_datasets=names_datasets,
                out_dir=out_dir,
            )
            generated["negative_control"] = list(result)
        except Exception:
            logger.exception("Error generating negative control plots")

    return generated
