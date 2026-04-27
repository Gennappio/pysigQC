#!/bin/bash
# Installation script for the sigQC R reference stack (Linux / macOS).
#
# Creates (or updates) a conda environment from environment-r.yml. Every R,
# CRAN, and Bioconductor dependency is pulled from conda-forge/bioconda, so
# no system libraries (libmagick++-dev, libgmp-dev, ...) and no BiocManager
# round-trips are required.
#
# Usage:
#   ./install_r.sh                # create/update env named "sigqc-r"
#   ./install_r.sh my-env-name    # use a custom env name
#   ./install_r.sh --update       # force `conda env update --prune` if env exists

set -e

cd "$(dirname "$0")"

ENV_NAME="sigqc-r"
FORCE_UPDATE=0
for arg in "$@"; do
    case "$arg" in
        --update) FORCE_UPDATE=1 ;;
        --help|-h)
            sed -n '2,12p' "$0"
            exit 0 ;;
        -*) echo "Unknown flag: $arg" >&2; exit 2 ;;
        *)  ENV_NAME="$arg" ;;
    esac
done

YML_FILE="environment-r.yml"

echo "======================================"
echo "sigQC R Dependencies Installation"
echo "  env: $ENV_NAME    spec: $YML_FILE"
echo "======================================"
echo ""

if [ ! -f "$YML_FILE" ]; then
    echo "ERROR: $YML_FILE not found in $(pwd)" >&2
    exit 1
fi

# ---- Step 1: Locate conda / mamba ----
echo "[STEP 1] Locating conda..."
if command -v mamba &> /dev/null; then
    CONDA_BIN="mamba"
elif command -v conda &> /dev/null; then
    CONDA_BIN="conda"
else
    echo "ERROR: neither 'conda' nor 'mamba' is on PATH." >&2
    echo "  Install Miniconda from https://docs.conda.io/en/latest/miniconda.html" >&2
    exit 1
fi
echo "  Using: $CONDA_BIN ($($CONDA_BIN --version | head -1))"
echo ""

# ---- Step 2: Create or update the conda env ----
# Detect whether the env already exists by name OR by full prefix.
ENV_EXISTS=0
if conda env list 2>/dev/null | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
    ENV_EXISTS=1
fi

if [ "$ENV_EXISTS" -eq 0 ]; then
    echo "[STEP 2] Creating conda env '$ENV_NAME' from $YML_FILE..."
    "$CONDA_BIN" env create -f "$YML_FILE" -n "$ENV_NAME"
elif [ "$FORCE_UPDATE" -eq 1 ]; then
    echo "[STEP 2] Updating conda env '$ENV_NAME' from $YML_FILE (--prune)..."
    "$CONDA_BIN" env update -f "$YML_FILE" -n "$ENV_NAME" --prune
else
    echo "[STEP 2] Conda env '$ENV_NAME' already exists - skipping (use --update to refresh)."
fi
echo ""

# Resolve the env's Rscript so we can run inside it without `conda activate`.
ENV_PREFIX="$($CONDA_BIN env list 2>/dev/null | awk -v n="$ENV_NAME" '$1==n {print $NF}')"
if [ -z "$ENV_PREFIX" ] || [ ! -x "$ENV_PREFIX/bin/Rscript" ]; then
    echo "ERROR: cannot locate Rscript in env '$ENV_NAME' (expected at $ENV_PREFIX/bin/Rscript)" >&2
    exit 1
fi
ENV_RSCRIPT="$ENV_PREFIX/bin/Rscript"
echo "Env prefix : $ENV_PREFIX"
echo "Rscript    : $ENV_RSCRIPT"
"$ENV_RSCRIPT" --version 2>&1 | head -1
echo ""

# ---- Step 3: Verify required packages load ----
echo "[STEP 3] Verifying R packages..."
"$ENV_RSCRIPT" - << 'EOF'
required <- c(
    "biclust", "fmsb", "moments", "mclust", "circlize", "gplots",
    "gridGraphics", "ComplexHeatmap", "GSVA", "biomaRt", "testthat",
    "flexclust", "additivityTests", "AnnotationDbi", "annotate"
)
optional <- c("RankProd", "devtools")

cat("\nRequired packages:\n")
missing <- c()
for (pkg in required) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat("  OK ", pkg, "-", as.character(packageVersion(pkg)), "\n")
    } else {
        cat("  -- ", pkg, "- MISSING\n")
        missing <- c(missing, pkg)
    }
}

cat("\nOptional packages:\n")
for (pkg in optional) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat("  OK ", pkg, "-", as.character(packageVersion(pkg)), "\n")
    } else {
        cat("  .. ", pkg, "- not installed (optional)\n")
    }
}

if (length(missing) > 0) {
    cat("\nMissing required packages:", paste(missing, collapse = ", "), "\n")
    quit(status = 1)
} else {
    cat("\nAll required packages available.\n")
}
EOF

echo ""
echo "======================================"
echo "Installation complete!"
echo "======================================"
echo ""
echo "Use the env without activating it by pointing RSCRIPT at:"
echo "  $ENV_RSCRIPT"
echo ""
echo "Next steps:"
echo "  1. Install sigQC package (env activated):  R CMD INSTALL sigQC-master/"
echo "  2. Run the cross-validation smoke test:"
echo "       RSCRIPT=\"$ENV_RSCRIPT\" \\"
echo "       SIGQC_USER_LIB=\"$ENV_PREFIX/lib/R/library\" \\"
echo "       python TEST_stud/test3.py"
echo ""
