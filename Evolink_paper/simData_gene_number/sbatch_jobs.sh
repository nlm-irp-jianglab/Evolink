#!/bin/bash
#SBATCH --partition largemem
#SBATCH --job-name simulate
#SBATCH --cpus-per-task=72
#SBATCH --time=5-00:00:00
source /data/$USER/conda/etc/profile.d/conda.sh
conda activate bio-env
cd $PWD/simData_gene_number
time Rscript --vanilla simdata_gene_number.R $1 $2 $3 $4
