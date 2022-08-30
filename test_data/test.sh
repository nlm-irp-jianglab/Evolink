# In this folder:
# test 1: no plot mode=gsed_test (default)
python ../Evolink.py -g gene.tsv -t trait.tsv -n tree.nwk -o test1 -f

# test 2: no plot mode=isolation_forest
python ../Evolink.py -g gene.tsv -t trait.tsv -n tree.nwk -o test2 -f -m isolation_forest

# test 3: no plot mode=z_score
python ../Evolink.py -g gene.tsv -t trait.tsv -n tree.nwk -o test3 -f -m z_score

# test 4: no plot mode=cutoff
python ../Evolink.py -g gene.tsv -t trait.tsv -n tree.nwk -o test4 -f -m cutoff

# test 5: no plot mode=permutation
python ../Evolink.py -g gene.tsv -t trait.tsv -n tree.nwk -o test5 -f -m permutation

# test 6: plot with default options
python ../Evolink.py -g gene.tsv -t trait.tsv -n tree.nwk -o test6 -v -f

# test 7: plot with top 1 pos and 1 neg genotypes
python ../Evolink.py -g gene.tsv -t trait.tsv -n tree.nwk -o test7 -v -N 1,1 -f

# test 8: no plot default mode with copy number gene matrix
python ../Evolink.py -g geneCN.tsv -c -t trait.tsv -n tree.nwk -o test8 -f

# clean folders
rm -rf test1 test2 test3 test4 test5 test6 test7 test8

# get example output
python ../Evolink.py -g gene.tsv -t trait.tsv -n tree.nwk -o test_output -f -v
