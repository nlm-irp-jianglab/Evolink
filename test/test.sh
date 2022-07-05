# In this folder:
# test 1: no plot no permutation
python ../Evolink.py -g gene.tsv -t trait.tsv -n tree.nwk -o test1 -f

# test 2: no plot with permutation
python ../Evolink.py -g gene.tsv -t trait.tsv -n tree.nwk -o test2 -s 100 -@ 8 -f

# test 3: plot but no permutation
python ../Evolink.py -g gene.tsv -t trait.tsv -n tree.nwk -o test3 -v -N 6,7 -f

# clean folders
rm -rf test1 test2 test3