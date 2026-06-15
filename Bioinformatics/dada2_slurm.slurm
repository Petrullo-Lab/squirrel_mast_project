#!/usr/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=dada2_16s
#SBATCH --output=outputr%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --account=laurenpetrullo
#SBATCH --mail-type=ALL
#SBATCH --mail-user=schiro@arizona.edu
#SBATCH --partition=standard  
#SBATCH --time=20:00:00
#SBATCH --ntasks=24
#SBATCH --nodes=1
#SBATCH --mem=120gb

# --------------------------------------------------
# Load modules here
# --------------------------------------------------
module load contrib
module load micromamba
eval "$(micromamba shell hook --shell bash)"  # Initialize micromamba for bash
micromamba activate /contrib/laurenpetrullo/dada2

# --------------------------------------------------
# Execute commands here
# --------------------------------------------------
mkdir   /groups/laurenpetrullo/data_sept_2025/fastq/dada2_outputs/
mkdir   /groups/laurenpetrullo/data_sept_2025/fastq/dada2_outputs/16s
Rscript /groups/laurenpetrullo/HPC_scripts/16s_dada2.R "/groups/laurenpetrullo/data_sept_2025/fastq/16S" \
 "/groups/laurenpetrullo/data_sept_2025/fastq/dada2_outputs/16s"