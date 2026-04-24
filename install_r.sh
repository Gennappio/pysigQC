#!/bin/bash
# Installation script for R and sigQC dependencies (Linux / macOS)
#
# This script:
#   1. Locates R — checks PATH, then common conda environments
#   2. Installs all CRAN dependencies
#   3. Installs biclust from archive (current CRAN version may have issues)
#   4a. Installs GSVA via conda (bioconda) when a conda env is active;
#       otherwise falls back to BiocManager (needs libmagick++-dev)
#   4b. Installs remaining Bioconductor packages via BiocManager
#       (ComplexHeatmap, biomaRt, AnnotationDbi, annotate, RankProd)
#   5. Installs testthat for running R tests

set -e

cd "$(dirname "$0")"

echo "======================================"
echo "sigQC R Dependencies Installation"
echo "======================================"
echo ""

# ---- Step 1: Locate R ----
echo "[STEP 1] Locating R..."

# If R is already on PATH, use it directly.
# Otherwise, search common conda/mamba environment locations.
if ! command -v R &> /dev/null; then
    echo "  R not found in PATH — searching conda environments..."

    R_SEARCH_PATHS=(
        "$HOME/miniconda3/envs/R_env/bin"
        "$HOME/miniconda3/envs/r_env/bin"
        "$HOME/anaconda3/envs/R_env/bin"
        "$HOME/anaconda3/envs/r_env/bin"
        "$HOME/mambaforge/envs/R_env/bin"
        "$HOME/mambaforge/envs/r_env/bin"
        "/opt/conda/envs/R_env/bin"
    )

    for dir in "${R_SEARCH_PATHS[@]}"; do
        if [ -x "$dir/R" ]; then
            export PATH="$dir:$PATH"
            echo "  Found R at: $dir/R"
            break
        fi
    done
fi

# Final check
if ! command -v R &> /dev/null; then
    echo "ERROR: R is not installed or not in PATH."
    echo "  Options:"
    echo "    a) Activate the conda environment first:  conda activate R_env"
    echo "    b) Install R from https://cran.r-project.org/"
    exit 1
fi

echo ""
echo "R version:"
R --version | head -1
echo ""

# ---- Step 1b: Install required system libraries ----
# Several R packages compile C/C++ extensions and need dev headers.
MISSING_SYS_LIBS=()
for pkg in libssl-dev libcurl4-openssl-dev libgmp-dev libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev libtiff5-dev; do
    if ! dpkg-query -W -f='${Status}' "$pkg" 2>/dev/null | grep -q "install ok installed"; then
        MISSING_SYS_LIBS+=("$pkg")
    fi
done

if [ ${#MISSING_SYS_LIBS[@]} -gt 0 ]; then
    echo "[STEP 1b] Missing system libraries: ${MISSING_SYS_LIBS[*]}"
    echo ""
    echo "  ERROR: The following system libraries must be installed before R packages"
    echo "  can compile. Please run the command below and then re-run this script:"
    echo ""
    echo "    sudo apt-get install -y ${MISSING_SYS_LIBS[*]}"
    echo ""
    exit 1
else
    echo "[STEP 1b] All required system libraries already present."
fi
echo ""

# ---- Ensure a writable R user library ----
# When the system library is not writable (e.g. conda R or system R without sudo),
# redirect package installation to a user-local directory.
R_USER_LIB="${HOME}/R/libs"
mkdir -p "$R_USER_LIB"
export R_LIBS_USER="$R_USER_LIB"
echo "R user library: $R_USER_LIB"
echo ""

# ---- Step 2: Install CRAN packages ----
echo "[STEP 2] Installing CRAN packages..."
Rscript - << 'EOF'
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
options(repos = c(CRAN = "https://cran.r-project.org/"))

# Core dependencies from DESCRIPTION
core_packages <- c(
    "BiocManager",
    "fmsb",          # radar charts
    "moments",       # skewness
    "mclust",        # Gaussian mixture models
    "circlize",      # circular plots (ComplexHeatmap dep)
    "gplots",        # heatmaps
    "gridGraphics",  # grid graphics capture
    "devtools"       # package development
)

# biclust dependencies (install separately before biclust)
biclust_deps <- c(
    "flexclust",
    "additivityTests",
    "tidyr",
    "ggplot2"
)

# Testing
test_packages <- c("testthat")

all_packages <- c(core_packages, biclust_deps, test_packages)

user_lib <- Sys.getenv("R_LIBS_USER")
cat("Installing", length(all_packages), "CRAN packages...\n")
for (pkg in all_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat("  Installing:", pkg, "\n")
        install.packages(pkg, lib = user_lib, quiet = TRUE)
    } else {
        cat("  Already installed:", pkg, "\n")
    }
}
cat("CRAN packages done.\n")
EOF

# ---- Step 3: Install biclust from archive ----
# The current CRAN version (2.0.4) may have issues; use archived 2.0.3.1
echo ""
echo "[STEP 3] Installing biclust (archived version 2.0.3.1)..."
Rscript - << 'EOF'
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
if (!requireNamespace("biclust", quietly = TRUE)) {
    cat("Installing biclust 2.0.3.1 from CRAN archive...\n")
    # Try archived version first (more stable)
    tryCatch({
        install.packages(
            "https://cran.r-project.org/src/contrib/Archive/biclust/biclust_2.0.3.1.tar.gz",
            repos = NULL,
            type = "source"
        )
        cat("biclust 2.0.3.1 installed from archive.\n")
    }, error = function(e) {
        cat("Archive install failed, trying current CRAN version...\n")
        install.packages("biclust", repos = "https://cran.r-project.org/")
    })
} else {
    cat("biclust already installed:", as.character(packageVersion("biclust")), "\n")
}
EOF

# ---- Step 4a: Install GSVA via conda (bioconda) ----
# GSVA >= 2.x depends on SpatialExperiment -> magick -> libmagick++-dev
# (system). Installing via BiocManager therefore requires sudo on the host.
# When a conda env is active, we pull bioconductor-gsva (plus imagemagick,
# magick, SpatialExperiment, ...) from the bioconda channel instead — the
# conda package bundles its own ImageMagick so no root is needed.
echo ""
echo "[STEP 4a] Installing GSVA via conda (bioconda)..."
if [ -n "$CONDA_PREFIX" ] && command -v conda &> /dev/null; then
    R_MAJMIN=$(R --version | head -1 | sed -E 's/.*version ([0-9]+)\.([0-9]+).*/\1\2/')
    # bioconda tags BioC 3.22 as r45 (R 4.5). For R 4.4 we would pick r44, etc.
    case "$R_MAJMIN" in
        45) GSVA_SPEC="bioconductor-gsva=2.4.4=r45h01b2380_0" ;;
        44) GSVA_SPEC="bioconductor-gsva=2.0.0" ;;
        43) GSVA_SPEC="bioconductor-gsva=1.50.0" ;;
        *)  GSVA_SPEC="bioconductor-gsva" ;;
    esac
    echo "  Conda env active: $CONDA_PREFIX"
    echo "  Installing: $GSVA_SPEC"
    conda install -y -n "$(basename "$CONDA_PREFIX")" \
        -c bioconda -c conda-forge "$GSVA_SPEC"
else
    echo "  No active conda env — will try BiocManager as fallback (may fail on"
    echo "  systems without libmagick++-dev). For a root-free install, activate"
    echo "  your R conda env first, e.g.:  conda activate R_env"
    export R_BIOC_VERSION="${R_BIOC_VERSION:-3.22}"
    Rscript - << 'EOF'
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cran.r-project.org/")
options(BiocManager.check_repositories = FALSE)
if (!requireNamespace("GSVA", quietly = TRUE)) {
    BiocManager::install("GSVA", ask = FALSE, update = FALSE)
} else {
    cat("  Already installed: GSVA\n")
}
EOF
fi

# ---- Step 4b: Install remaining Bioconductor packages via BiocManager ----
# These do not have the magick/SpatialExperiment chain and install cleanly
# from BiocManager on any platform.
echo ""
echo "[STEP 4b] Installing remaining Bioconductor packages..."
export R_BIOC_VERSION="${R_BIOC_VERSION:-3.22}"
Rscript - << 'EOF'
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cran.r-project.org/")

bioc_packages <- c(
    "annotate",       # used by ComplexHeatmap/biomaRt stacks
    "AnnotationDbi",  # biomaRt dependency
    "ComplexHeatmap",
    "biomaRt",        # Ensembl/HGNC gene ID mapping (used by TEST_stud/test3.R)
    "RankProd"        # optional; needs libgmp-dev system library
)

options(BiocManager.check_repositories = FALSE)
for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat("  Installing:", pkg, "\n")
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
        cat("  Already installed:", pkg, "\n")
    }
}
cat("Bioconductor packages done.\n")
EOF

# ---- Step 5: Verify installation ----
echo ""
echo "[STEP 5] Verifying all packages..."
Rscript - << 'EOF'
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
required <- c(
    "biclust", "fmsb", "moments", "mclust", "circlize", "gplots",
    "gridGraphics", "ComplexHeatmap", "GSVA", "biomaRt", "testthat"
)
optional <- c("RankProd")

cat("\nRequired packages:\n")
missing <- c()
for (pkg in required) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat("  ✓", pkg, "-", as.character(packageVersion(pkg)), "\n")
    } else {
        cat("  ✗", pkg, "- MISSING\n")
        missing <- c(missing, pkg)
    }
}

cat("\nOptional packages:\n")
for (pkg in optional) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat("  ✓", pkg, "-", as.character(packageVersion(pkg)), "\n")
    } else {
        cat("  ○", pkg, "- not installed (optional)\n")
    }
}

if (length(missing) > 0) {
    cat("\n⚠️  Missing required packages:", paste(missing, collapse = ", "), "\n")
    quit(status = 1)
} else {
    cat("\n✅ All required packages installed successfully!\n")
}
EOF

echo ""
echo "======================================"
echo "Installation complete!"
echo "======================================"
echo ""
echo "Next steps:"
echo "  1. Install sigQC package:  R CMD INSTALL sigQC-master/"
echo "  2. Generate test fixtures: cd sigQC-master && Rscript tests/fixtures/fixture_generator.R"
echo "  3. Run R tests:            cd sigQC-master && Rscript tests/testthat.R"
echo ""
