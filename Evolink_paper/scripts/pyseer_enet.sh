#!/bin/bash

########################################################
# Wrap Pyseer using elastic net model in a bash script #
########################################################

source /data/$USER/conda/etc/profile.d/conda.sh
conda activate pyseer_env

tree_path=$1
gene_path=$2
trait_path=$3
outdir=$4
alpha=${5:-0.0069}

mkdir -p ${outdir}

# gene matrix
sed "s/#OTU ID/gene/g" ${gene_path} > ${outdir}/gene.tsv
gene_path=${outdir}/gene.tsv

# trait 
awk -F"\t" 'NR==1{print "samples\tbinary"}NR>1{print $0}' ${trait_path} > ${outdir}/tmp
mv ${outdir}/tmp ${outdir}/trait.tsv
trait_path=${outdir}/trait.tsv

# distance matrix
python softwares/PySeer/phylogeny_distance.py ${tree_path} > ${outdir}/phylogeny_dist.tsv
dist_path=${outdir}/phylogeny_dist.tsv

# enet
pyseer --pres ${gene_path} --phenotypes ${trait_path} --wg enet --distances ${dist_path} --save-vars ${outdir}/vars_alpha${alpha} --cpu 32 --alpha ${alpha} > ${outdir}/selected_alpha${alpha}.txt
cat ${outdir}/selected_alpha${alpha}.txt > ${outdir}/output_enet_patterns_alpha${alpha}.txt
threshold=`python softwares/PySeer/count_patterns.py --alpha 0.05 --cores 32 --memory 1000 --temp /tmp ${outdir}/output_enet_patterns_alpha${alpha}.txt |awk 'NR==2'|cut -f2`
echo ${threshold}
awk -F"\t" -v thresh=${threshold} '$4<thresh{print $1}' ${outdir}/output_enet_patterns_alpha${alpha}.txt > ${outdir}/Pyseer_enet_alpha${alpha}_sig_genes.list
