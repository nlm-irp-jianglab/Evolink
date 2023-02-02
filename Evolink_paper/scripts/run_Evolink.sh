#!/bin/bash

tree_file=$1
gene_file=$2
trait_file=$3
outdir=$4
mode=${5:-isolation_forest}

mkdir -p ${outdir}

source /data/$USER/conda/etc/profile.d/conda.sh
conda activate bio-env

python softwares/Evolink.py -n ${tree_file} -g ${gene_file} -t ${trait_file} -o ${outdir} -m ${mode} -f
grep -w sig ${outdir}/result.tsv| cut -f1 > ${outdir}/Evolink_sig_genes.list
