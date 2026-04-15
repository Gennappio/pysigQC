#!/bin/bash
# Installation script for R and sigQC dependencies on macOS
#
# This script:
#   1. Installs R 4.5.3 (if .pkg file present) OR uses existing R
#   2. Installs all CRAN dependencies
#   3. Installs biclust from archive (current CRAN version may have issues)
#   4. Installs Bioconductor packages (ComplexHeatmap, GSVA, RankProd)
#   5. Installs testthat for running R tests

set -e

cd "$(dirname "$0")"

echo "======================================"
echo "sigQC R Dependencies Installation"
echo "======================================"
echo ""

# ---- Step 1: Install R if .pkg file present ----
if [ -f "R-4.5.3-arm64.pkg" ]; then
    echo "[STEP 1] Installing R 4.5.3 for Apple Silicon (arm64)..."
    echo "This requires administrator (sudo) access."
    sudo installer -pkg R-4.5.3-arm64.pkg -target /
    echo "R installation complete!"
else
    echo "[STEP 1] R-4.5.3-arm64.pkg not found, using existing R installation."
fi

# Verify R is available
if ! command -v R &> /dev/null; then
    echo "ERROR: R is not installed or not in PATH."
    echo "Install R from https://cran.r-project.org/ or provide R-4.5.3-arm64.pkg"
    exit 1
fi

echo ""
echo "R version:"
R --version | head -1
echo ""

# ---- Step 2: Install CRAN packages ----
echo "[STEP 2] Installing CRAN packages..."
Rscript << 'EOF'
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

cat("Installing", length(all_packages), "CRAN packages...\n")
for (pkg in all_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat("  Installing:", pkg, "\n")
        install.packages(pkg, quiet = TRUE)
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
Rscript << 'EOF'
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

# ---- Step 4: Install Bioconductor packages ----
echo ""
echo "[STEP 4] Installing Bioconductor packages..."
Rscript << 'EOF'
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cran.r-project.org/")

bioc_packages <- c("ComplexHeatmap", "GSVA", "RankProd")

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
Rscript << 'EOF'
required <- c(
    "biclust", "fmsb", "moments", "mclust", "circlize", "gplots",
    "gridGraphics", "ComplexHeatmap", "GSVA", "testthat"
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
