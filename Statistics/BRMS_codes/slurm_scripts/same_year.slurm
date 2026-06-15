#!/usr/bin/bash
#SBATCH --job-name=brms_same_year
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --account=laurenpetrullo
#SBATCH --mail-type=ALL
#SBATCH --mail-user=schiro@arizona.edu
#SBATCH --partition=standard  
#SBATCH --time=50:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120gb

# -------------------------------
# Environment setup
# -------------------------------
module load R

# Prevent over-threading inside BLAS/OpenMP libraries
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Move to working directory
cd /xdisk/laurenpetrullo/schiro/dyadic_models/new_subselection

# -------------------------------
# Run R script
# -------------------------------
Rscript BRMS_same_year_diff_excluded.R