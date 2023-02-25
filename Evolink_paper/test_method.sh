#!/bin/bash

wd_dir=$1 # a folder containing gene_file, trait file, and tree file
method=$2

scrdir=scripts

gene_file=${simulation_dir}/gene.tsv
trait_file=${simulation_dir}/trait.tsv
tree_file=${simulation_dir}/tree.nwk

cd ${wd_dir}

case $method in

Tet_corr)
    # tet_corr
    echo "Tetrachoric correlation"
    Tet_corr=${scrdir}/run_tetcorr.R
    mkdir -p Tet_corr_out
    time Rscript --vanilla ${Tet_corr} ${gene_file} ${trait_file} Tet_corr_out &> Tet_corr.out
;;

FG)
    # FG
    echo "ForwardGenomics"
    ForwardGenomics_GLS=${scrdir}/run_ForwardGenomics.sh
    time /bin/bash ${ForwardGenomics_GLS} ${tree_file} ${gene_file} ${trait_file} FG_out &> FG.out
;;

Phylolm)
    # Phylolm
    echo "Phylolm"
    Phylolm=${scrdir}/run_Phylolm.sh
    time /bin/bash ${Phylolm} ${tree_file} ${gene_file} ${trait_file} Phylolm_out &> Phylolm.out
;;

Pyseer_fem)
    # Pyseer
    # fem
    echo "Pyseer fem"
    Pyseer_fem=${scrdir}/pyseer_fem.sh
    time /bin/bash ${Pyseer_fem} ${tree_file} ${gene_file} ${trait_file} Pyseer_fem_out &> Pyseer_fem.out
;;

Pyseer_lmm)
    # lmm
    echo "Pyseer lmm"
    Pyseer_lmm=${scrdir}/pyseer_lmm.sh
    time /bin/bash ${Pyseer_lmm} ${tree_file} ${gene_file} ${trait_file} Pyseer_lmm_out &> Pyseer_lmm.out
;;

Pyseer_enet_alpha0)
    # ridge reg
    echo "Pyseer enet alpha0"
    Pyseer_enet=${scrdir}/pyseer_enet.sh
    time /bin/bash ${Pyseer_enet} ${tree_file} ${gene_file} ${trait_file} Pyseer_enet_alpha0_out 0 &> Pyseer_enet_alpha0.out
;;

Pyseer_enet_alpha0.5)
    # middle
    echo "Pyseer enet alpha0.5"
    Pyseer_enet=${scrdir}/pyseer_enet.sh
    time /bin/bash ${Pyseer_enet} ${tree_file} ${gene_file} ${trait_file} Pyseer_enet_alpha0.5_out 0.5 &> Pyseer_enet_alpha0.5.out
;;

Pyseer_enet_alpha1)
    # lasso reg
    echo "Pyseer enet alpha1"
    Pyseer_enet=${scrdir}/pyseer_enet.sh
    time /bin/bash ${Pyseer_enet} ${tree_file} ${gene_file} ${trait_file} Pyseer_enet_alpha1_out 1 &> Pyseer_enet_alpha1.out
;;

BUGWAS)
    # BUGWAS
    echo "BUGWAS"
    BUGWAS=${scrdir}/run_BUGWAS.sh
    time /bin/bash ${BUGWAS} ${tree_file} ${gene_file} ${trait_file} BUGWAS_out &> BUGWAS.out
;;

Hogwash_phyc)
    # Hogwash
    echo "Hogwash"
    Hogwash=${scrdir}/run_hogwash_phyc.sh
    time /bin/bash ${Hogwash} ${tree_file} ${gene_file} ${trait_file} Hogwash_phyc_out &> Hogwash_phyc.out
;;

Hogwash_synchronous)
    # Hogwash
    echo "Hogwash"
    Hogwash=${scrdir}/run_hogwash_synchronous.sh
    time /bin/bash ${Hogwash} ${tree_file} ${gene_file} ${trait_file} Hogwash_synchronous_out &> Hogwash_synchronous.out
;;

treeWAS)
    # treeWAS
    echo "treeWAS"
    treeWAS=${scrdir}/run_treeWAS.R
    time Rscript ${treeWAS} ${tree_file} ${gene_file} ${trait_file} treeWAS_out &> treeWAS.out
;;

Evolink)
    # Evolink
    echo "Evolink"
    Evolink=${scrdir}/run_Evolink.sh
    time /bin/bash ${Evolink} ${tree_file} ${gene_file} ${trait_file} Evolink_out &> Evolink.out
;;

esac
