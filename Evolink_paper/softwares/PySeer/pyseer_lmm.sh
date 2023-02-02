#!/bin/bash

source /data/$USER/conda/etc/profile.d/conda.sh
conda activate pyseer_env

tree_path=$1
gene_path=$2
trait_path=$3
outdir=$4

mkdir -p ${outdir}

# gene matrix
sed "s/#OTU ID/gene/g" ${gene_path} > ${outdir}/gene.tsv
gene_path=${outdir}/gene.tsv

# trait 
awk -F"\t" 'NR==1{print "samples\tbinary"}NR>1{print $0}' ${trait_path} > tmp
mv tmp ${outdir}/trait.tsv
trait_path=${outdir}/trait.tsv

# distance matrix
python /data/yangy34/softwares/PySeer/phylogeny_distance.py --lmm ${tree_path} > ${outdir}/phylogeny_K.tsv
dist_path=${outdir}/phylogeny_K.tsv

# lmm, mixed effects model
pyseer --lmm --phenotypes ${trait_path} --pres ${gene_path} --similarity ${dist_path} --output-patterns ${outdir}/pres_patterns.txt --cpu 32 > ${outdir}/output_lmm.txt
cat ${outdir}/output_lmm.txt |sort -g -k4,4 -t$'\t' > ${outdir}/output_lmm_patterns.txt
threshold=`python /data/yangy34/softwares/PySeer/count_patterns.py --alpha 0.05 --cores 32 --memory 1000 --temp /tmp ${outdir}/output_lmm_patterns.txt|awk 'NR==2'|cut -f2`
echo ${threshold}
awk -F"\t" -v thresh=${threshold} '$4<thresh{print $1}' ${outdir}/output_lmm_patterns.txt > ${outdir}/sig_gene.list
