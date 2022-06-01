#!/bin/bash
#SBATCH --partition norm
#SBATCH --job-name v2_Evolink
#SBATCH --cpus-per-task=56
#SBATCH --mem=240g
#SBATCH --time=8-00:00:00
source /data/$USER/conda/etc/profile.d/conda.sh
conda activate bio-env
cd /data/yangy34/projects/Evolink/test_flagella
python /data/yangy34/projects/Evolink/Evolink_v2.py -g gene.tsv -t trait.tsv -n tree.nwk -o flagella_genes_v2.tsv -p 1000
