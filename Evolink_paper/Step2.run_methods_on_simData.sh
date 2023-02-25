dirs=(simData_phenotype_phylo_overdispersion simData_genotype_number simData_phenotype_prevalence simData_species_number)
methods=(Tet_corr FG Phylolm Pyseer_fem Pyseer_lmm Pyseer_enet_alpha0 Pyseer_enet_alpha0.5 Pyseer_enet_alpha1 BUGWAS Hogwash_phyc Hogwash_synchronous treeWAS Evolink)
mkdir -p logs
for dir in ${dirs[*]}
do
    for wkdir in $(ls -d ${dir}/*/)
    do  
        wd=$(realpath ${wkdir})
        for method in ${methods[*]}
        do
            var=${method}.$(basename ${wkdir})
            echo $var
            # This step is to run each method on a HPC node by submitting it to Slurm
            sbatch --error=logs/${var}.err --output=logs/${var}.out run_job.sh ${wd} ${method}
            # echo "sbatch --error=logs/${var}.err --output=logs/${var}.out run_job.sh ${wd} ${method}"
        done
    done
done
