#!/bin/bash
#SBATCH --partition largemem
#SBATCH --job-name simulate
#SBATCH --cpus-per-task=72
#SBATCH --mem=800g
#SBATCH --time=5-00:00:00
source /data/$USER/conda/etc/profile.d/conda.sh
conda activate bio-env
cd $PWD/simData_species_number
time Rscript --vanilla simdata_species_number.R $1 $2 $3 $4
