#!/bin/bash
#SBATCH --job-name=brms_asv
#SBATCH --account=laurenpetrullo
#SBATCH --mail-type=ALL
#SBATCH --mail-user=schiro@arizona.edu
#SBATCH --partition=standard  
#SBATCH --array=1-50
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120gb
#SBATCH --time=32:00:00
#SBATCH --output=logs/genera_%A_%a.out
#SBATCH --error=logs/genera_%A_%a.err

# Go to working directory
cd /xdisk/laurenpetrullo/schiro/dyadic_models/ord_beta_reg
echo "directory set"

# Make sure output directories exist
mkdir -p output_models
mkdir -p logs

# Load R module
module load R
echo "R module loaded at $(date)"
which Rscript
Rscript --version

# Show which array index is running
echo "Running array index: $SLURM_ARRAY_TASK_ID"

# Run the R script, passing array index
Rscript code_genera_new.R $SLURM_ARRAY_TASK_ID