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
awk -F"\t" 'NR==1{print "samples\tbinary"}NR>1{print $0}' ${trait_path} > ${outdir}/tmp
mv ${outdir}/tmp ${outdir}/trait.tsv
trait_path=${outdir}/trait.tsv

# distance matrix
python softwares/PySeer/phylogeny_distance.py ${tree_path} > ${outdir}/phylogeny_dist.tsv
dist_path=${outdir}/phylogeny_dist.tsv

# fixed effects model
pyseer --phenotypes ${trait_path} --pres ${gene_path} --distances ${dist_path} --cpu 32 > ${outdir}/output_fem.txt
cat ${outdir}/output_fem.txt |sort -g -k4,4 -t$'\t' > ${outdir}/output_fem_patterns.txt # sort by lrt-pvalue column
threshold=`python softwares/PySeer/count_patterns.py --alpha 0.05 --cores 32 --memory 1000 --temp /tmp ${outdir}/output_fem_patterns.txt|awk 'NR==2'|cut -f2`
echo ${threshold}
awk -F"\t" -v thresh=${threshold} '$4<thresh{print $1}' ${outdir}/output_fem_patterns.txt > ${outdir}/Pyseer_fem_sig_genes.list
