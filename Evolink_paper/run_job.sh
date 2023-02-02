#!/bin/bash
#SBATCH --partition norm
#SBATCH --job-name test_on_simData
#SBATCH --cpus-per-task=32
#SBATCH --mem=120g
#SBATCH --time=8-00:00:00
#SBATCH --gres=lscratch:100
source /data/$USER/conda/etc/profile.d/conda.sh
conda activate /data/yangy34/conda/envs/bio-env

export LD_LIBRARY_PATH="/usr/lib64:$LD_LIBRARY_PATH"
export PATH="/usr/lib64:$PATH"

dir=$1
method=$2
/bin/bash test_method.sh ${dir} ${method}
