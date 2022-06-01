#!/bin/bash
#SBATCH --partition largemem
#SBATCH --job-name v1_Evolink
#SBATCH --cpus-per-task=144
#SBATCH --mem=350g
#SBATCH --time=3-00:00:00
source /data/$USER/conda/etc/profile.d/conda.sh
conda activate bio-env
cd /data/yangy34/projects/Evolink/test_flagella
export OMP_NUM_THREADS=4
python /data/yangy34/projects/Evolink/Evolink_v1.py -g gene.tsv -t trait.tsv -n tree.nwk -o flagella_genes_v1.tsv -p 1000 -@ 36
