#!/bin/bash
#SBATCH --job-name="sigQC_original"
#SBATCH --mem=16000
#SBATCH --account=abbruzzese
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=6
#SBATCH --output=out/%x_%j.out
#SBATCH --error=err/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=gennaro.abbruzzese@unibocconi.it

echo "start"
echo "Running original sigQC (unmodified) on BRCA paper example data"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start time: $(date)"

cd "$(dirname "$0")/../input_files"

Rscript run_original_with_timing.R

echo "End time: $(date)"
echo "done"
