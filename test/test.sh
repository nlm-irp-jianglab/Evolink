# In this folder:
# test 1: no plot no permutation
python ../Evolink.py -g gene.tsv -t trait.tsv -n tree.nwk -o test1 -f

# test 2: no plot with permutation
python ../Evolink.py -g gene.tsv -t trait.tsv -n tree.nwk -o test2 -s 1000 -@ 8 -f

# test 3: plot but no permutation
python ../Evolink.py -g gene.tsv -t trait.tsv -n tree.nwk -o test3 -v -N 6,7 -f

# test 4: plot and permutation
python ../Evolink.py -g gene.tsv -t trait.tsv -n tree.nwk -o test4 -s 1000 -@ 8 -v -N 6,7 -f

# test 5: with given Evolink index threshold
python ../Evolink.py -g gene.tsv -t trait.tsv -n tree.nwk -o test5 -f -e 0.45 -v

# clean folders
rm -rf test1 test2 test3 test4 test5

# get example output
python ../Evolink.py -g gene.tsv -t trait.tsv -n tree.nwk -o output_dir -f -v
python ../Evolink.py -g gene.tsv -t trait.tsv -n tree.nwk -o output_perm_dir -f -s 1000 -@ 8 -v
