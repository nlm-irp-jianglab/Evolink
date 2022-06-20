#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os, tempfile, datetime
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from unifrac import faith_pd
from biom.util import biom_open
from biom.table import Table
from tqdm import tqdm

def Evolink_calculation(trait_matrix, gene_matrix, species_list, gene_list, tree_file):

    T1_index, T0_index, G1T1, G0T1, G1T0, G0T0 = generate_four_matrix(gene_matrix, trait_matrix)

    T1_species = [species_list[i] for i in T1_index]
    T0_species = [species_list[i] for i in T0_index]
    T1_path, T1_len = generate_tree(tree_file,T1_species)
    T0_path, T0_len = generate_tree(tree_file,T0_species)

    G1T1_obs = matrix_faith_pd(G1T1, T1_species, gene_list, T1_path)
    G0T1_obs = matrix_faith_pd(G0T1, T1_species, gene_list, T1_path)
    G1T0_obs = matrix_faith_pd(G1T0, T0_species, gene_list, T0_path)
    G0T0_obs = matrix_faith_pd(G0T0, T0_species, gene_list, T0_path)
    # delete temp tree files
    os.unlink(T1_path)
    os.unlink(T0_path)

    T1_obs = np.tile(round(T1_len, 5), len(gene_list))
    T0_obs = np.tile(round(T0_len, 5), len(gene_list))
    
    G1T1_T1 = G1T1_obs/T1_obs
    G0T1_T1 = G0T1_obs/T1_obs
    G1T0_T0 = G1T0_obs/T0_obs
    G0T0_T0 = G0T0_obs/T0_obs
    
    Evolink_index = (G1T1_T1 - G1T0_T0 + G0T0_T0 - G0T1_T1)*0.5
    Prevalence_index = (G1T1_T1 + G1T0_T0 - G0T1_T1 - G0T0_T0)*0.5

    result = np.vstack((Prevalence_index, Evolink_index)).T
    df = pd.DataFrame(result, columns=["Prevalence_index", "Evolink_index"], index=gene_list)
    return df

def matrix2hdf5(mat, row_names, col_names, outfile):
    mat = Table(mat, row_names, col_names)
    with biom_open(outfile, 'w') as f:
        mat.to_hdf5(f, "Evolink")

def matrix_faith_pd( matrix, species_name, gene_list, tree_path):
    fd, path = tempfile.mkstemp(suffix=".biom")
    matrix2hdf5(matrix, species_name, gene_list, path)
    obs = faith_pd(path,tree_path)
    os.close(fd)
    # delete temp biom files
    os.unlink(path)
    return obs

def generate_four_matrix(gene_matrix, trait_matrix):
    T1_index = np.where(trait_matrix == [1])[0]
    T0_index = np.where(trait_matrix == [0])[0]
    G1T1 = np.multiply(gene_matrix,np.transpose(trait_matrix)).T[T1_index]
    G0T1 = np.multiply(1-gene_matrix,np.transpose(trait_matrix)).T[T1_index]
    G1T0 = np.multiply(gene_matrix,np.transpose(1-trait_matrix)).T[T0_index]
    G0T0 = np.multiply(1-gene_matrix,np.transpose(1-trait_matrix)).T[T0_index]
    return (T1_index, T0_index, G1T1, G0T1, G1T0, G0T0)

def sum_br_len(tree_file):
    res = robjects.r('''
    library(ape)
    tree = read.tree("'''+tree_file+'''")
    sum(tree$edge.length)
    ''')
    br_len_sum = res[0]
    return br_len_sum

def generate_tree(tree_file, species_sub):
    fd, path = tempfile.mkstemp(suffix=".nwk")
    species = robjects.vectors.StrVector(species_sub)
    robjects.globalenv['species'] = species
    robjects.r('''
    library(ape)
    tree = read.tree("'''+tree_file+'''")
    pruned.tree<-drop.tip(tree, setdiff(tree$tip.label, species));
    write.tree(pruned.tree, "'''+path+'''")
    ''')
    br_len = sum_br_len(path)
    os.close(fd)
    return (path, br_len)

def sig_genes(df, p_threshold, e_threshold, alpha):
    with localconverter(robjects.default_converter + pandas2ri.converter):
        r_df = robjects.conversion.py2rpy(df)
    robjects.globalenv['r_df'] = r_df
    robjects.globalenv['p_threshold'] = p_threshold
    robjects.globalenv['e_threshold'] = e_threshold
    robjects.globalenv['alpha'] = alpha
    rcode = '''
    suppressPackageStartupMessages({
        library(fdrtool)
        suppressWarnings(library(tidyverse))
    })
    r_df = r_df %>% rownames_to_column(var = "orthoID")
    dt = r_df %>% filter(Prevalence_index >= -1*p_threshold)
    fdrres = fdrtool(dt$Evolink_index, statistic="normal", plot=F, verbose=F)
    dt$pvalue=fdrres$pval
    dt$qvalue=fdrres$qval
    dt = dt %>% mutate(significance=ifelse(pvalue < 0.05 & qvalue < alpha & abs(Evolink_index) >= e_threshold, "sig", NA))
    r_df = r_df %>% select(orthoID, Prevalence_index, Evolink_index) %>% 
                left_join(dt, by = c("orthoID", "Prevalence_index", "Evolink_index")) %>% 
                select(orthoID, Prevalence_index, Evolink_index, significance, pvalue, qvalue)
    '''
    r_df = robjects.r(rcode)
    with localconverter(robjects.default_converter + pandas2ri.converter):
        res_df = robjects.conversion.rpy2py(r_df)
    return res_df

def pipeline(args):

    gene_table = args.gene_table
    trait_table = args.trait_table
    tree_file = args.tree_file
    CN = args.CN
    p_threshold = args.p_threshold
    e_threshold = args.e_threshold
    alpha = args.alpha
    output = args.output

    if CN:
        cn = "-c "
    else:
        cn = ""

    ### Print command line ###
    print("Command line: python Evolink.py -g {0} -t {1} -n {2} {3}-f {4} -e {5} -a {6} -o {7}".format(gene_table, trait_table,
                                                                                                        tree_file, cn,
                                                                                                        p_threshold, e_threshold,
                                                                                                        alpha, output), flush=True)
    ### Read trait data ###
    print("[",datetime.datetime.now(),"]","Read trait data", flush=True)
    trait_df = pd.read_csv(trait_table, sep='\t',index_col=0)
    species_list = trait_df.index.values
    trait_matrix = trait_df.values

    ### Read gene data ###
    print("[",datetime.datetime.now(),"]","Read gene data", flush=True)
    gene_df = pd.read_csv(gene_table, sep='\t',index_col=0)
    gene_df = gene_df[species_list]
    if CN:
        gene_df[:] = np.where(gene_df > 0, 1, 0)
    gene_list = gene_df.index.values
    gene_species_list = gene_df.columns.values
    gene_matrix = gene_df.values

    # test if the species in gene_table is the same to those in trait list
    if set(gene_species_list) == set(species_list):
        pass
    else:
        raise ValueError("Dimmension of genotype file do not match that of trait file")

    ### Calculate Evolink index ###
    print("[",datetime.datetime.now(),"]","Calculate Evolink index", flush=True)
    df = Evolink_calculation(trait_matrix, gene_matrix, 
                             species_list, gene_list, tree_file)

    ### Get pvalues and qvalues ###
    print("[",datetime.datetime.now(),"]","Find significant genes", flush=True)
    res_df = sig_genes(df, p_threshold, e_threshold, alpha)

    ### Output ###
    print("[",datetime.datetime.now(),"]","Output result", flush=True)
    res_df.to_csv(output, sep="\t", index=False, na_rep='NA')

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Evolink is designed to find gene families associated with given trait with the help of phylogeny information.\n'+ \
    'Examples: [With binary gene table] python Evolink.py -g test/gene.tsv -t test/trait.tsv -n test/tree.nwk -o test_res.tsv; \n'+ \
    '[With non-binary gene table] python Evolink.py -g test/gene_CN.tsv -c -t test/trait.tsv -n test/tree.nwk -o test_res.tsv.')
    
    # Essential Input
    parser.add_argument('-g', '--genotype', help='Tab-delimited gene presence/absence or copy number table. Columns are gene families, while rows are tip names/species/genomes in the phylogenetic tree. If copy number table is provided, please use -c option so that it will be internally converted to binary table. Presence=1, Absence=0.', required=True, dest='gene_table', metavar='GENE_TABLE')
    parser.add_argument('-t', '--phenotype', help='Two-column (so far only one trait is allowed each time) tab-delimited trait presence/absence table. The first column is tip names and the second column is the presence/absence of this trait on the tips/species/genomes. Presence=1, Absence=0.',required=True, dest='trait_table', metavar='TRAIT_TABLE')
    parser.add_argument('-n', '--phylogeny', help='A phylogentic tree in newick format. The tip names should be the same in the gene table and trait table.',required=True, dest='tree_file', metavar='TREE')

    # Optional Input
    parser.add_argument('-c', '--copy_number', help='The given gene table stores numbers (e.g. gene copy numbers) instead of presence/absence binary values. [Default: True]', action='store_true', required=False, default=False, dest='CN')
    parser.add_argument('-f', '--p_threshold', help='Absolute Prevalence index cutoff to filter in genes for permutation tests [Range: 0-1; Default: 0.9]', required=False, default=0.9, type=float, dest='p_threshold', metavar='THRESHOLD')
    parser.add_argument('-e', '--e_threshold', help='Absolute Evolink index cutoff to filter in genes for permutation tests [Range: 0-1; Default: 0.1]', required=False, default=0.1, type=float, dest='e_threshold', metavar='THRESHOLD')
    parser.add_argument('-a', '--alpha', help='Adjusted p-value (or qvalue) cutoff [Range: 0-1; Default: 0.05]', required=False, default=0.05, type=float, dest='alpha', metavar='ALPHA')
    
    # Output
    parser.add_argument('-o', '--output', help='Output file', required=True, dest='output', metavar='OUTPUT')
    
    args = parser.parse_args()
    pipeline(args)
