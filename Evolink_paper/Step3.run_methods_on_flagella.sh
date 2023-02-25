# run on flagella full dataset
mkdir -p test_flagella_data
unzip DataS5_flagella_data.zip -d test_flagella_data/
methods=(Phylolm Pyseer_fem Pyseer_lmm Pyseer_enet_alpha0 Pyseer_enet_alpha0.5 Pyseer_enet_alpha1 BUGWAS treeWAS Evolink)
mkdir -p logs
for method in ${methods[*]}
do
    wd=test_flagella_data
    var=${method}.test_flagella_data
    echo $var
    # This step is to run method on a HPC node by submitting it to Slurm
    sbatch --error=logs/${var}.err --output=logs/${var}.out run_job.sh ${wd} ${method}
done

# run on flagella small dataset for ForwardGenomics and Hogwash
mkdir -p test_flagella_subdata
unzip ../flagella_data/example_input_subset.zip -d test_flagella_subdata/
methods=(FG Hogwash_phyc Hogwash_synchronous)
mkdir -p logs
for method in ${methods[*]}
do
    wd=test_flagella_subdata
    var=${method}.test_flagella_subdata
    echo $var
    # This step is to run method on a HPC node by submitting it to Slurm
    sbatch --error=logs/${var}.err --output=logs/${var}.out run_job.sh ${wd} ${method}
done
