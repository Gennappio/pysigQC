#!/usr/bin/env python
"""
pysigQC Command-Line Interface

Run gene signature quality control on genomic datasets from the command line.

Example usage:
    # Single dataset and signature file
    python main.py --datasets data/expr.csv --signatures sigs.csv --out-dir results/

    # Multiple datasets with glob pattern
    python main.py --datasets "data/dataset_*.csv" --signatures sigs.csv --out-dir results/

    # Disable plotting (compute only)
    python main.py --datasets data.csv --signatures sigs.csv --out-dir results/ --no-plots

    # Change input format (parquet, tsv, etc.)
    python main.py --datasets data.parquet --signatures sigs.csv --format parquet

    # Control parallelism
    python main.py --datasets data.csv --signatures sigs.csv --n-jobs 4
"""

from __future__ import annotations

import argparse
import glob
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


def load_dataset(path: Path, fmt: str) -> pd.DataFrame:
    """Load a single expression matrix (genes x samples)."""
    if fmt == "auto":
        suffix = path.suffix.lower()
        fmt = {".csv": "csv", ".tsv": "tsv", ".parquet": "parquet", ".h5": "hdf"}.get(suffix, "csv")

    loaders = {
        "csv": lambda p: pd.read_csv(p, index_col=0),
        "tsv": lambda p: pd.read_csv(p, index_col=0, sep="\t"),
        "parquet": lambda p: pd.read_parquet(p).set_index(pd.read_parquet(p).columns[0]),
        "hdf": lambda p: pd.read_hdf(p),
    }
    return loaders[fmt](path)


def load_signatures(path: Path, sig_format: str) -> dict[str, list[str]]:
    """Load gene signatures from file.

    Supports:
      - csv: long format with columns 'signature' and 'gene'
      - json: dict of signature_name -> [gene1, gene2, ...]
      - gmt: MSigDB GMT format (name<tab>description<tab>gene1<tab>gene2...)
    """
    if sig_format == "auto":
        suffix = path.suffix.lower()
        sig_format = {".csv": "csv", ".tsv": "csv", ".json": "json", ".gmt": "gmt"}.get(suffix, "csv")

    sigs: dict[str, list[str]] = {}

    if sig_format == "csv":
        df = pd.read_csv(path)
        for sig_name, group in df.groupby("signature"):
            sigs[str(sig_name)] = group["gene"].astype(str).tolist()

    elif sig_format == "json":
        with open(path) as f:
            data = json.load(f)
        sigs = {k: list(v) for k, v in data.items()}

    elif sig_format == "gmt":
        with open(path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    name, _desc, *genes = parts
                    sigs[name] = genes

    return sigs


def _save_tables(
    result: dict,
    names_sigs: list[str],
    names_datasets: list[str],
    out_dir: Path,
    logger: logging.Logger,
) -> None:
    """Save intermediate CSV tables matching R's output folder structure.

    R produces these subdirectories:
      radarchart_table/        — the final summary table
      mean_sd_tables/          — per-gene mean and SD for signature genes
      expression_tables/       — per-gene NA proportion and expression proportion
      autocorrelation_matrices/ — gene-gene Spearman correlation matrices
      metrics_tables/          — PCA loadings and variance explained
      standardisation_tables/  — raw vs z-transformed median scores per sample
    """
    # --- radarchart_table ---
    radar_dir = out_dir / "radarchart_table"
    radar_dir.mkdir(exist_ok=True)
    radar_df = result["radar_result"]["output_table"]
    radar_df.to_csv(radar_dir / "radarchart_table.txt", sep="\t")
    # Also save at top level for convenience
    radar_df.to_csv(out_dir / "radar_output_table.csv")
    logger.info("Saved radarchart_table/")

    # --- mean_sd_tables ---
    var_result = result.get("var_result", {})
    msd = var_result.get("mean_sd_tables", {})
    if msd:
        d = out_dir / "mean_sd_tables"
        d.mkdir(exist_ok=True)
        for sig in names_sigs:
            for ds in names_datasets:
                tbl = msd.get(sig, {}).get(ds)
                if tbl is not None:
                    tbl.to_csv(d / f"mean_sd_table_{sig}_{ds}.txt", sep="\t")
        logger.info("Saved mean_sd_tables/")

    # --- expression_tables ---
    expr_result = result.get("expr_result", {})
    na_props = expr_result.get("na_proportions", {})
    expr_props = expr_result.get("expr_proportions", {})
    if na_props:
        d = out_dir / "expression_tables"
        d.mkdir(exist_ok=True)
        for sig in names_sigs:
            for ds in names_datasets:
                na_s = na_props.get(sig, {}).get(ds)
                ex_s = expr_props.get(sig, {}).get(ds)
                if na_s is not None and ex_s is not None:
                    tbl = pd.DataFrame({"prop_NA": na_s, "prop_expressed": ex_s})
                    tbl.to_csv(d / f"expression_table_{sig}_{ds}.txt", sep="\t")
        logger.info("Saved expression_tables/")

    # --- autocorrelation_matrices ---
    compact_result = result.get("compact_result", {})
    autocor = compact_result.get("autocor_matrices", {})
    if autocor:
        d = out_dir / "autocorrelation_matrices"
        d.mkdir(exist_ok=True)
        for sig in names_sigs:
            for ds in names_datasets:
                mat = autocor.get(sig, {}).get(ds)
                if mat is not None:
                    mat.to_csv(d / f"autocorrelation_matrix_{sig}_{ds}.txt", sep="\t")
        logger.info("Saved autocorrelation_matrices/")

    # --- metrics_tables (PCA loadings + variance) ---
    metrics_result = result.get("metrics_result", {})
    scores = metrics_result.get("scores", {})
    pca_results = metrics_result.get("pca_results", {})
    if scores:
        d = out_dir / "metrics_tables"
        d.mkdir(exist_ok=True)
        for sig in names_sigs:
            for ds in names_datasets:
                sc = scores.get(sig, {}).get(ds, {})
                props = sc.get("props_of_variances")
                if props is not None and len(props) > 0:
                    pca_var = pd.DataFrame(
                        {"PC": range(1, len(props) + 1), "prop_variance": props}
                    )
                    pca_var.to_csv(
                        d / f"pca_vs_var_{sig}_{ds}.txt", sep="\t", index=False
                    )
                pca_obj = pca_results.get(sig, {}).get(ds, {}).get("pca_obj")
                if pca_obj is not None and hasattr(pca_obj, "components_"):
                    import numpy as np
                    loadings = pd.DataFrame(
                        pca_obj.components_.T,
                        columns=[f"PC{i+1}" for i in range(pca_obj.n_components_)],
                    )
                    cols = sc.get("common_score_cols")
                    if cols is not None and len(cols) == loadings.shape[0]:
                        loadings.index = cols
                    loadings.to_csv(
                        d / f"pca_loadings_{sig}_{ds}.txt", sep="\t"
                    )
        logger.info("Saved metrics_tables/")

    # --- standardisation_tables ---
    stan_result = result.get("stan_result", {})
    med = stan_result.get("med_scores", {})
    z_med = stan_result.get("z_transf_scores", {})
    if med:
        d = out_dir / "standardisation_tables"
        d.mkdir(exist_ok=True)
        for sig in names_sigs:
            for ds in names_datasets:
                m = med.get(sig, {}).get(ds)
                z = z_med.get(sig, {}).get(ds)
                if m is not None and z is not None:
                    import numpy as np
                    tbl = pd.DataFrame({"median_score": np.asarray(m),
                                        "z_median_score": np.asarray(z)})
                    tbl.to_csv(
                        d / f"standardisation_table_{sig}_{ds}.txt",
                        sep="\t", index=False,
                    )
        logger.info("Saved standardisation_tables/")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="pysigQC: Quality control for gene signatures across genomic datasets.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    # Required inputs
    parser.add_argument(
        "--datasets", "-d", required=True, nargs="+",
        help="Path(s) or glob pattern(s) to expression matrices (genes x samples). "
             "First column/index = gene IDs, remaining columns = samples.",
    )
    parser.add_argument(
        "--signatures", "-s", required=True,
        help="Path to gene signatures file (CSV with 'signature','gene' columns, JSON dict, or GMT).",
    )
    parser.add_argument(
        "--out-dir", "-o", default="sigqc_output",
        help="Output directory for results and plots (default: sigqc_output).",
    )

    # Format options
    parser.add_argument(
        "--format", "-f", default="auto", choices=["auto", "csv", "tsv", "parquet", "hdf"],
        help="Input format for expression matrices (default: auto-detect from extension).",
    )
    parser.add_argument(
        "--sig-format", default="auto", choices=["auto", "csv", "json", "gmt"],
        help="Signature file format (default: auto-detect).",
    )

    # Pipeline control
    parser.add_argument(
        "--no-plots", action="store_true",
        help="Skip plot generation (compute metrics only).",
    )
    parser.add_argument(
        "--n-jobs", "-j", type=int, default=-1,
        help="Number of parallel jobs (-1 = all cores, default: -1).",
    )
    parser.add_argument(
        "--save-results", action="store_true", default=True,
        help="Save intermediate results as JSON/CSV (default: True).",
    )

    # Logging
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose output.")
    parser.add_argument("--quiet", "-q", action="store_true", help="Suppress non-error output.")

    args = parser.parse_args(argv)
    setup_logging(args.verbose and not args.quiet)
    logger = logging.getLogger("pysigqc")

    if args.quiet:
        logging.getLogger().setLevel(logging.ERROR)

    # -------------------------------------------------------------------------
    # Load datasets
    # -------------------------------------------------------------------------
    dataset_paths: list[Path] = []
    for pattern in args.datasets:
        matches = glob.glob(pattern)
        if not matches:
            logger.warning("No files matched pattern: %s", pattern)
        dataset_paths.extend(Path(m) for m in sorted(matches))

    if not dataset_paths:
        logger.error("No dataset files found. Check --datasets paths/patterns.")
        return 1

    logger.info("Loading %d dataset(s)...", len(dataset_paths))
    mRNA_expr_matrix: dict[str, pd.DataFrame] = {}
    for p in dataset_paths:
        name = p.stem  # Use filename without extension as dataset name
        logger.debug("  Loading %s as '%s'", p, name)
        mRNA_expr_matrix[name] = load_dataset(p, args.format)

    names_datasets = list(mRNA_expr_matrix.keys())
    logger.info("Datasets: %s", names_datasets)

    # -------------------------------------------------------------------------
    # Load signatures
    # -------------------------------------------------------------------------
    sig_path = Path(args.signatures)
    if not sig_path.exists():
        logger.error("Signature file not found: %s", sig_path)
        return 1

    logger.info("Loading signatures from %s...", sig_path)
    gene_sigs_list = load_signatures(sig_path, args.sig_format)
    names_sigs = list(gene_sigs_list.keys())
    logger.info("Signatures (%d): %s", len(names_sigs), names_sigs)

    # -------------------------------------------------------------------------
    # Run pipeline
    # -------------------------------------------------------------------------
    from pysigqc_joblib.pipeline import run_pipeline

    logger.info("Running QC pipeline (n_jobs=%s)...", args.n_jobs)
    result = run_pipeline(
        gene_sigs_list=gene_sigs_list,
        names_sigs=names_sigs,
        mRNA_expr_matrix=mRNA_expr_matrix,
        names_datasets=names_datasets,
        n_jobs=args.n_jobs,
        verbose=args.verbose,
    )
    logger.info("Pipeline completed in %.2f seconds.", result["elapsed_seconds"])

    # -------------------------------------------------------------------------
    # Save results
    # -------------------------------------------------------------------------
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.save_results:
        _save_tables(result, names_sigs, names_datasets, out_dir, logger)

    # -------------------------------------------------------------------------
    # Generate plots
    # -------------------------------------------------------------------------
    if not args.no_plots:
        try:
            from pysigqc.plots import plot_all

            logger.info("Generating plots...")
            plot_paths = plot_all(result, out_dir=out_dir)
            total_plots = sum(len(v) for v in plot_paths.values())
            logger.info("Generated %d plot files in %s", total_plots, out_dir)
        except ImportError:
            logger.warning("Plotting dependencies not installed. Run: pip install pysigqc[plot]")
        except Exception as e:
            logger.exception("Error generating plots: %s", e)

    logger.info("Done. Results saved to %s", out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
