# run on gram-staining full dataset
mkdir -p test_gram_staining_data
unzip DataS6_gram_staining_data.zip -d test_gram_staining_data/
methods=(Phylolm Pyseer_fem Pyseer_lmm Pyseer_enet_alpha0 Pyseer_enet_alpha0.5 Pyseer_enet_alpha1 BUGWAS Evolink)
mkdir -p logs
for method in ${methods[*]}
do
    wd=test_gram_staining_data
    var=${method}.test_gram_staining_data
    echo $var
    # This step is to run method on a HPC node by submitting it to Slurm
    sbatch --error=logs/${var}.err --output=logs/${var}.out run_job.sh ${wd} ${method}
done

# run on gram staining small dataset for treeWAS
mkdir -p test_gram_staining_subdata
unzip ../gram_staining_data/example_input_subset.zip -d test_gram_staining_subdata/
methods=(treeWAS)
mkdir -p logs
for method in ${methods[*]}
do
    wd=test_gram_staining_subdata
    var=${method}.test_gram_staining_subdata
    echo $var
    # This step is to run method on a HPC node by submitting it to Slurm
    sbatch --error=logs/${var}.err --output=logs/${var}.out run_job.sh ${wd} ${method}
done
