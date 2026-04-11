#!/bin/bash
#SBATCH --job-name="sigQC_py_tests"
#SBATCH --mem=8000
#SBATCH --account=abbruzzese
#SBATCH --partition=long_gpunew
#SBATCH --cpus-per-task=4
#SBATCH --output=out/%x_%j.out
#SBATCH --error=err/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=gennaro.abbruzzese@unibocconi.it

set -e
echo "Job ID: $SLURM_JOB_ID  |  Node: $SLURMD_NODENAME  |  Start: $(date)"

SCRIPT_DIR="$(realpath "${SLURM_SUBMIT_DIR}")"
PROJECT_DIR="$(realpath "${SCRIPT_DIR}/../..")"
ENV_NAME="sigqc_py_env"

# ---- Initialize conda ----
source "$(conda info --base)/etc/profile.d/conda.sh"

# ---- Create or reuse conda environment ----
if ! conda env list | grep -q "^${ENV_NAME} "; then
  echo "[INSTALL] Creating Python conda environment..."
  conda create -y -n "${ENV_NAME}" python=3.11
  echo "[INSTALL] Conda environment created."
else
  echo "[INSTALL] Conda environment '${ENV_NAME}' already exists."
fi

conda activate "${ENV_NAME}"

# ---- Install pysigqc in dev mode ----
echo "[INSTALL] Installing pysigqc[dev]..."
cd "${PROJECT_DIR}"
pip install -e ".[dev]" --quiet

echo "[INSTALL] Done."

# ---- Verify fixtures exist ----
FIXTURES_DIR="${PROJECT_DIR}/sigQC-master/tests/fixtures"
if [ ! -f "${FIXTURES_DIR}/fixture_dataset_A.csv" ]; then
  echo "[ERROR] Fixtures not found in ${FIXTURES_DIR}."
  echo "        Run run_r_tests.slurm first to generate fixtures and reference outputs."
  exit 1
fi

# ---- Run Python unit tests ----
echo ""
echo "========== PYTHON UNIT TESTS (pytest) =========="
echo "Working directory: $(pwd)"
echo ""

echo "[TEST] Running unit tests..."
pytest tests/ -v --tb=short -x 2>&1

echo ""
echo "========== PYTHON TESTS COMPLETE =========="
echo "End time: $(date)"
