#!/bin/bash

tree_file=$1
gene_file=$2
trait_file=$3
outdir=$4

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
mkdir -p ${outdir}

source /data/$USER/conda/etc/profile.d/conda.sh
conda activate bio-env

awk -F"\t" 'BEGIN{OFS="\t"}NR==1{$1="gene"; print $0}NR>1{print $0}' ${gene_file} > ${outdir}/gene4hogwash.tsv
Rscript ${SCRIPTPATH}/run_hogwash_phyc.R ${tree_file} ${outdir}/gene4hogwash.tsv ${trait_file} ${outdir} &
wait
# rm *.pdf
mv *.rda ${outdir}
Rscript ${SCRIPTPATH}/hogwash_phyc_sig_gene.R ${outdir}
