#!/bin/bash

source /data/$USER/conda/etc/profile.d/conda.sh
conda activate bio-env

ml LAPACK/3.8.0/gcc-4.8.5
export LD_LIBRARY_PATH="/usr/lib64:$LD_LIBRARY_PATH"
export PATH="/usr/lib64:$PATH"

tree_file=$1
gene_file=$2
trait_file=$3
outprefix=$4

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
mkdir -p ${outprefix}

sed "s/\b0\b/C/g;s/\b1\b/A/g" ${gene_file} | awk -F"\t" 'BEGIN{OFS="\t"}NR==1{$1="ps"; print $0}NR>1{print $0}' > ${outprefix}/gene4bugwas.tmp
# remove invariant sites in gene file
Rscript ${SCRIPTPATH}/BUGWAS_modGene.R ${outprefix}/gene4bugwas.tmp ${outprefix}/gene4bugwas.tsv
awk 'BEGIN{OFS="\t"}NR==1{print "ID\tpheno"}NR>1{print $0}' ${trait_file} > ${outprefix}/trait4bugwas.tsv

Rscript --vanilla ${SCRIPTPATH}/run_BUGWAS.R ${tree_file} ${outprefix}/gene4bugwas.tsv ${outprefix}/trait4bugwas.tsv bugwas_output &
wait
mv bugwas_output* ${outprefix}
rm -rf output/
Rscript ${SCRIPTPATH}/BUGWAS_sig_gene.R ${outprefix}
# rm ${outprefix}/bugwas_output*
