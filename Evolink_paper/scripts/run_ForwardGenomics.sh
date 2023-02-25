#!/bin/bash

#################################################
# Wrap ForwardGenomics running in a bash script #
#################################################

source /data/$USER/conda/etc/profile.d/conda.sh
conda activate bio-env

tree_file=$1
gene_file=$2
trait_file=$3
outdir=$4

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

mkdir -p ${outdir}

Rscript ${SCRIPTPATH}/FG_modTree.R ${tree_file} ${outdir}/tree4FG.nwk

awk -F"\t" 'NR==1{print "species pheno"}NR>1{print $1" "$2}' ${trait_file} > ${outdir}/trait4FG.txt

sed "s/#OTU ID/species/g;s/\t/ /g" ${gene_file} > ${outdir}/gene4FG.txt
sed -i "s/\(^[0-9]\{1\}[0-9A-Z].*\)/X\1/g" ${outdir}/gene4FG.txt
mv ${outdir}/gene4FG.txt ${outdir}/globalpid.txt
tail -n +2 ${outdir}/globalpid.txt | cut -f1 -d " " > ${outdir}/id.list

export PATH=$PATH:$PWD/softwares/phast/usr/bin

# only GLS
softwares/ForwardGenomics/forwardGenomics.R --tree=${outdir}/tree4FG.nwk --elementIDs=${outdir}/id.list --listPheno=${outdir}/trait4FG.txt --globalPid=${outdir}/globalpid.txt --method=GLS --outFile=${outdir}/output.txt
Rscript ${SCRIPTPATH}/FG_sig_gene.R ${outdir}/output.txt ${outdir}