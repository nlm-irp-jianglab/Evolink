# run on gram-staining full dataset
mkdir -p test_gram_staining_data
unzip DataS6_gram_staining_data.zip -d test_gram_staining_data/
methods=(Phylolm Pyseer_fem Pyseer_lmm Pyseer_enet_alpha0 Pyseer_enet_alpha0.5 Pyseer_enet_alpha1 BUGWAS treeWAS Evolink)
for method in ${methods[*]}
do
    time test_method.sh test_gram_staining_data ${method} &> ${method}.log
done
