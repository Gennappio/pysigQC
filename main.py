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
                # Combined scoring metrics per sample
                med = sc.get("med_scores")
                mean = sc.get("mean_scores")
                pca1 = sc.get("pca1_scores")
                col_names = sc.get("common_score_cols")
                if med is not None and mean is not None and col_names:
                    import numpy as np
                    metrics_df = pd.DataFrame(
                        {"Mean_Scores": mean, "Median_Scores": med},
                        index=col_names,
                    )
                    if pca1 is not None and len(pca1) == len(col_names):
                        metrics_df["PCA1_Scores"] = pca1
                    # Include enrichment scores if available
                    for es_key, es_col in [
                        ("ssgsea_scores", "ssGSEA"),
                        ("plage_scores", "PLAGE"),
                        ("gsva_scores", "GSVA"),
                    ]:
                        es = sc.get(es_key)
                        if es is not None and len(es) == len(col_names):
                            metrics_df[es_col] = es
                    metrics_df.to_csv(
                        d / f"metrics_table_{sig}_{ds}.txt", sep="\t"
                    )
        logger.info("Saved metrics_tables/")

    # --- mixture_model_out.txt ---
    mixture_models = metrics_result.get("mixture_models", {})
    if mixture_models:
        lines = []
        for sig in names_sigs:
            for ds in names_datasets:
                mm = mixture_models.get(sig, {}).get(ds, {})
                for score_name in ("median", "mean", "pca1"):
                    entry = mm.get(score_name)
                    if isinstance(entry, dict) and entry.get("best_model"):
                        best_k = entry["best_k"]
                        if best_k == 1:
                            lines.append(
                                f"{sig} {ds}, {score_name.title()} score: "
                                f"There is one component in the mixture model."
                            )
                        else:
                            lines.append(
                                f"{sig} {ds}, {score_name.title()} score: "
                                f"There are {best_k} components in the mixture model."
                            )
        if lines:
            (out_dir / "mixture_model_out.txt").write_text("\n".join(lines) + "\n")
            logger.info("Saved mixture_model_out.txt")

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

    # --- negative/permutation control tables ---
    neg_result = result.get("negative_result")
    if neg_result:
        for control_type in ("negative_controls", "permutation_controls"):
            controls = neg_result.get(control_type, {})
            subdir = "negative_control" if "negative" in control_type else "permutation_control"
            for ds in names_datasets:
                for sig in names_sigs:
                    data = controls.get(ds, {}).get(sig)
                    if data is None:
                        continue
                    d = out_dir / subdir / ds / sig
                    d.mkdir(parents=True, exist_ok=True)
                    summary = data.get("summary")
                    if summary is not None:
                        summary.to_csv(
                            d / f"{subdir.rstrip('_control')}_controls_summary_table.txt",
                            sep="\t",
                        )
                    metrics_table = data.get("metrics_table")
                    if metrics_table is not None:
                        sigqc_dir = d / "sigQC" / "radarchart_table"
                        sigqc_dir.mkdir(parents=True, exist_ok=True)
                        metrics_table.to_csv(
                            sigqc_dir / "radarchart_table.txt", sep="\t",
                        )
        logger.info("Saved negative/permutation control tables")


def _validate_inputs(
    gene_sigs_list: dict[str, list[str]],
    names_sigs: list[str],
    mRNA_expr_matrix: dict[str, pd.DataFrame],
    names_datasets: list[str],
    logger: logging.Logger,
) -> tuple[dict, list, dict, list]:
    """Validate and filter inputs, matching R's make_all_plots() validation.

    Removes signatures with <2 genes (after intersection) and datasets with
    <2 samples. Warns about low gene overlap. Exits if nothing remains.
    """
    from pysigqc.utils import gene_intersection

    # Remove datasets with <2 samples
    valid_datasets = {}
    for ds in names_datasets:
        n_samples = mRNA_expr_matrix[ds].shape[1]
        if n_samples < 2:
            logger.warning("Removing dataset '%s': only %d sample(s) (need >=2)", ds, n_samples)
        else:
            valid_datasets[ds] = mRNA_expr_matrix[ds]
    if not valid_datasets:
        logger.error("No datasets with >=2 samples. Aborting.")
        raise SystemExit(1)
    mRNA_expr_matrix = valid_datasets
    names_datasets = list(valid_datasets.keys())

    # Remove signatures with <2 genes in ANY dataset
    valid_sigs = {}
    for sig in names_sigs:
        genes = gene_sigs_list[sig]
        keep = True
        for ds in names_datasets:
            inter = gene_intersection(genes, mRNA_expr_matrix[ds])
            if len(inter) < 2:
                logger.warning(
                    "Removing signature '%s': only %d gene(s) found in dataset '%s' (need >=2)",
                    sig, len(inter), ds,
                )
                keep = False
                break
        if keep:
            valid_sigs[sig] = genes
    if not valid_sigs:
        logger.error("No signatures with >=2 genes in all datasets. Aborting.")
        raise SystemExit(1)
    gene_sigs_list = valid_sigs
    names_sigs = list(valid_sigs.keys())

    # Warn about low gene overlap
    for sig in names_sigs:
        genes = gene_sigs_list[sig]
        for ds in names_datasets:
            inter = gene_intersection(genes, mRNA_expr_matrix[ds])
            pct = len(inter) / len(genes) * 100
            if pct < 50:
                logger.warning(
                    "Signature '%s' in dataset '%s': only %.0f%% of genes found (%d/%d). "
                    "Check gene ID format.",
                    sig, ds, pct, len(inter), len(genes),
                )

    return gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets


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
    parser.add_argument(
        "--negative-control", action="store_true",
        help="Run negative and permutation controls (slow, default: off).",
    )
    parser.add_argument(
        "--num-resampling", type=int, default=50,
        help="Number of resamplings for negative control (default: 50).",
    )
    parser.add_argument(
        "--thresholds", type=float, nargs="+", default=None,
        help="Expression thresholds per dataset (default: median of each dataset).",
    )
    parser.add_argument(
        "--covariates", type=str, default=None,
        help="Path to covariates JSON file for heatmap annotations. "
             "Format: {dataset: {annotations: {sample: {var1: val, ...}}, colors: {var1: {val: color}}}}",
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
    # Validate inputs (matching R's make_all_plots() checks)
    # -------------------------------------------------------------------------
    gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets = _validate_inputs(
        gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets, logger
    )

    # -------------------------------------------------------------------------
    # Run pipeline
    # -------------------------------------------------------------------------
    from pysigqc_joblib.pipeline import run_pipeline

    # Validate thresholds if provided
    thresholds = args.thresholds
    if thresholds is not None:
        if len(thresholds) != len(names_datasets):
            logger.error(
                "Number of thresholds (%d) must match number of datasets (%d)",
                len(thresholds), len(names_datasets)
            )
            return 1
        logger.info("Using custom thresholds: %s", dict(zip(names_datasets, thresholds)))

    # Load covariates if provided
    covariates = None
    if args.covariates:
        cov_path = Path(args.covariates)
        if not cov_path.exists():
            logger.error("Covariates file not found: %s", cov_path)
            return 1
        with open(cov_path) as f:
            cov_data = json.load(f)
        # Convert annotations to DataFrames
        covariates = {}
        for ds, ds_cov in cov_data.items():
            covariates[ds] = {
                "annotations": pd.DataFrame(ds_cov.get("annotations", {})).T,
                "colors": ds_cov.get("colors", {}),
            }
        logger.info("Loaded covariates for datasets: %s", list(covariates.keys()))

    logger.info("Running QC pipeline (n_jobs=%s)...", args.n_jobs)
    result = run_pipeline(
        gene_sigs_list=gene_sigs_list,
        names_sigs=names_sigs,
        mRNA_expr_matrix=mRNA_expr_matrix,
        names_datasets=names_datasets,
        n_jobs=args.n_jobs,
        verbose=args.verbose,
        do_negative_control=args.negative_control,
        num_resampling=args.num_resampling,
        thresholds=thresholds,
        covariates=covariates,
    )
    logger.info("Pipeline completed in %.2f seconds.", result["elapsed_seconds"])

    # -------------------------------------------------------------------------
    # Save results
    # -------------------------------------------------------------------------
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Set up log file (matching R's log.log)
    log_path = out_dir / "log.log"
    file_handler = logging.FileHandler(log_path, mode="w")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s")
    )
    logging.getLogger().addHandler(file_handler)
    logger.info("pysigQC log started")

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
