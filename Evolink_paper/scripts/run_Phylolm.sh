#!/bin/bash

source /data/$USER/conda/etc/profile.d/conda.sh
conda activate bio-env

tree_file=$1
gene_file=$2
trait_file=$3
outdir=$4

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

mkdir -p ${outdir}

awk -F"\t" 'BEGIN{OFS="\t"}NR==1{$1="gene"; print $0}NR>1{print $0}' ${gene_file} > ${outdir}/gene4phylolm.tsv
Rscript --vanilla ${SCRIPTPATH}/run_phylolm.R -n ${tree_file} -g ${outdir}/gene4phylolm.tsv -t ${trait_file} -o ${outdir}/output.txt &
wait
Rscript ${SCRIPTPATH}/phylolm_sig_gene.R ${outdir}/output.txt ${outdir}