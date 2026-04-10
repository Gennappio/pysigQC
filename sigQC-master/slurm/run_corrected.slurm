#!/bin/bash
#SBATCH --job-name="sigQC_corrected"
#SBATCH --mem=16000
#SBATCH --account=abbruzzese
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=6
#SBATCH --output=out/%x_%j.out
#SBATCH --error=err/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=gennaro.abbruzzese@unibocconi.it

echo "start"
echo "Running corrected sigQC (8 bug fixes applied) on BRCA paper example data"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start time: $(date)"

cd /home/Abbruzzese/pysigQC/sigQC-master/input_files

Rscript run_corrected_with_timing.R

echo "End time: $(date)"
echo "done"
