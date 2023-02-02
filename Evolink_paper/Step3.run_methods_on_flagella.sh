# run on flagella full dataset
mkdir -p test_flagella_data
unzip DataS5_flagella_data.zip -d test_flagella_data/
methods=(Phylolm Pyseer_fem Pyseer_lmm Pyseer_enet_alpha0 Pyseer_enet_alpha0.5 Pyseer_enet_alpha1 BUGWAS treeWAS Evolink)
for method in ${methods[*]}
do
    time test_method.sh test_flagella_data ${method} &> ${method}.log
done

# run on flagella small dataset for ForwardGenomics and Hogwash
mkdir -p test_flagella_subdata
unzip ../flagella_data/example_input_subset.zip -d test_flagella_subdata/
methods=(FG Hogwash_phyc Hogwash_synchronous)
for method in ${methods[*]}
do
    time test_method.sh test_flagella_subdata ${method} &> ${method}.log
done
