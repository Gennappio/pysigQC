#!/usr/bin/env python
"""
pysigQC Plotting CLI — Generate plots from pre-computed results.

Use this when you've already run computations (e.g., with --no-plots)
and want to generate visualizations separately.

Example usage:
    # Plot from a saved radar_output_table.csv
    python main_plot.py --results results/radar_output_table.csv --out-dir plots/

    # Plot from a full pipeline result JSON
    python main_plot.py --results results/pipeline_result.json --out-dir plots/

    # Generate only specific plot types
    python main_plot.py --results results/ --out-dir plots/ --only radar var
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import pandas as pd


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )


def load_results(path: Path) -> dict:
    """Load pipeline results from file or directory.
    
    Supports:
      - JSON file with full pipeline result dict
      - Directory containing radar_output_table.csv and module CSVs
      - Single CSV file (radar_output_table.csv)
    """
    if path.is_file():
        if path.suffix == ".json":
            with open(path) as f:
                return json.load(f)
        elif path.suffix == ".csv":
            # Just radar table - create minimal result dict
            df = pd.read_csv(path, index_col=0)
            return {"radar_result": {"output_table": df}}
    
    elif path.is_dir():
        result = {}
        radar_csv = path / "radar_output_table.csv"
        if radar_csv.exists():
            result["radar_result"] = {
                "output_table": pd.read_csv(radar_csv, index_col=0)
            }
        return result
    
    raise ValueError(f"Cannot load results from: {path}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate plots from pre-computed pysigQC results.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--results", "-r", required=True,
        help="Path to results file (JSON) or directory containing CSVs.",
    )
    parser.add_argument(
        "--out-dir", "-o", default="plots",
        help="Output directory for plots (default: plots).",
    )
    parser.add_argument(
        "--only", nargs="+",
        choices=["var", "expr", "compactness", "stan", "struct", "metrics", "radar", "negative_control"],
        help="Generate only specific plot types.",
    )
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output.")

    args = parser.parse_args(argv)
    setup_logging(args.verbose)
    logger = logging.getLogger("pysigqc.plot")

    # Load results
    results_path = Path(args.results)
    if not results_path.exists():
        logger.error("Results path not found: %s", results_path)
        return 1

    logger.info("Loading results from %s...", results_path)
    try:
        result = load_results(results_path)
    except Exception as e:
        logger.error("Failed to load results: %s", e)
        return 1

    # Generate plots
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        from pysigqc.plots import plot_all
        
        logger.info("Generating plots in %s...", out_dir)
        plot_paths = plot_all(result, out_dir=out_dir)
        
        total = sum(len(v) for v in plot_paths.values())
        logger.info("Generated %d plot files:", total)
        for module, paths in plot_paths.items():
            for p in paths:
                logger.info("  %s: %s", module, p.name)
                
    except ImportError:
        logger.error("Plotting dependencies not installed. Run: pip install pysigqc[plot]")
        return 1
    except Exception as e:
        logger.exception("Error generating plots: %s", e)
        return 1

    logger.info("Done. Plots saved to %s", out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
