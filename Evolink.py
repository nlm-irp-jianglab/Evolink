#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os, tempfile, shutil

def Evolink_calculation(species_list,gene_species_list,gene_list,T1_index,T0_index, gene_matrix,trait_matrix,tree_file):

    T1_species = [species_list[i] for i in T1_index]
    T0_species = [species_list[i] for i in T0_index]

    G1T1, G0T1, G1T0, G0T0 = compute_four_matrix(gene_species_list, species_list,T1_species,T0_species, gene_matrix, trait_matrix)

    T1_path, T1_len = generate_tree(tree_file,T1_species)
    T0_path, T0_len = generate_tree(tree_file,T0_species)

    G1T1_obs = matrix_faith_pd(G1T1,T1_species,gene_list,T1_path)
    G0T1_obs = matrix_faith_pd(G0T1,T1_species,gene_list,T1_path)
    G1T0_obs = matrix_faith_pd(G1T0,T0_species,gene_list,T0_path)
    G0T0_obs = matrix_faith_pd(G0T0,T0_species,gene_list,T0_path)
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
    
    result = np.vstack((G1T1_obs, G0T1_obs, T1_obs, G1T0_obs, G0T0_obs, T0_obs, Prevalence_index, Evolink_index)).T
    df = pd.DataFrame(result, columns=["G1T1", "G0T1", "T1", "G1T0", "G0T0", "T0", "Prevalence_index", "Evolink_index"], index=gene_list)

    return df 

def matrix_faith_pd( matrix, species_name, gene_list,tree_path):
    from unifrac import faith_pd
    fd, path = tempfile.mkstemp(suffix=".biom")
    matrix2hdf5(matrix, species_name, gene_list, path)
    obs = faith_pd(path,tree_path)
    os.close(fd)
    os.unlink(path)
    return obs

def compute_four_matrix(gene_species_list, species_list,T1_species,T0_species, gene_matrix, trait_matrix):
    gene_perm_matrix = create_permuation_matrix(gene_species_list,species_list)
    T1_perm_matrix = create_permuation_matrix(species_list,T1_species)
    T0_perm_matrix = create_permuation_matrix(species_list,T0_species)

    gene_matrix_perm = gene_matrix @ gene_perm_matrix 

    G1T1 =( gene_matrix_perm     * trait_matrix       @ T1_perm_matrix).T 
    G0T1 =( (1-gene_matrix_perm) * trait_matrix       @ T1_perm_matrix).T 
    G1T0 =( gene_matrix_perm     * (1-trait_matrix)   @ T0_perm_matrix).T 
    G0T0 =( (1-gene_matrix_perm) * (1-trait_matrix)   @ T0_perm_matrix).T 

    return (G1T1, G0T1, G1T0, G0T0)

def generate_tree(tree_file,species_sub):
    from ete3 import Tree
    in_tree = Tree(tree_file,format=1)
    in_tree.prune(species_sub, preserve_branch_length = True)
    br_len = sum_br_len(in_tree)
    fd, path = tempfile.mkstemp(suffix=".nwk")
    in_tree.write(format=1,outfile=path,format_root_node=True)
    return (path, br_len)

def sum_br_len(tree):
    return sum([node.dist for node in tree.traverse() if not node.is_root()])


def matrix2hdf5(mat, row_names, col_names, outfile):
    from biom.util import biom_open
    from biom.table import Table

    mat = Table(mat, row_names, col_names)
    with biom_open(outfile, 'w') as f:
        mat.to_hdf5(f, "Evolink")

def create_permuation_matrix(list1,list2):
    matrix= np.zeros((len(list1),len(list2)))
    for element in list2:
        idx1 = list1.index(element)
        idx2 = list2.index(element)
        matrix[idx1,idx2] = 1
    return matrix

def pipeline(args):


    gene_table = args.gene_table
    trait_table = args.trait_table
    tree_file = args.tree_file
    CN = args.CN
    output = args.output

    trait_df = pd.read_csv(trait_table, sep='\t',index_col=0)
    # species_list is the one that can be permutated
    species_list = trait_df.index.values.tolist()
    trait_matrix = np.ravel(trait_df.values)
    T1_index = np.where(trait_matrix == 1)[0]
    T0_index = np.where(trait_matrix == 0)[0]

    gene_df = pd.read_csv(gene_table, sep='\t',index_col=0)
    if CN:
        gene_df.values = np.where(gene_df.values > 0, 1, 0)
    gene_list = gene_df.index.values
    gene_species_list = gene_df.columns.values.tolist()
    gene_matrix = gene_df.values

    if set(gene_species_list) == set(species_list):
        pass
    else:
        raise ValueError("dimmension of genotype file do not match that of trait file")

    # function : gene_species_list, species_list, T1_index, T0_index, gene_matrix, trait_matrix, tree_file,gene_list
    permutation_test = 10
    count = 0
    while count < permutation_test:
        species_list_sim = np.random.permutation(species_list).tolist()
        df_sim = Evolink_calculation(species_list_sim,gene_species_list,gene_list,T1_index,T0_index, gene_matrix,trait_matrix,tree_file)
        df_sim["Evolink_index"].to_csv(output+str(count),sep="\t")
        count += 1


    # species_list_sim = np.random.permutation(species_list)
    # Each permutation change the species_list order 
    
    df = Evolink_calculation(species_list,gene_species_list,gene_list,T1_index,T0_index, gene_matrix,trait_matrix,tree_file)

    df.to_csv(output, sep="\t")

    # clear intermediate files
    #shutil.rmtree(dirpath)

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Evolink is designed to find gene families associated with given trait with the help of phylogeny information.\n'+ \
    'Examples: [With binary gene table] python Evolink.py -g test/gene.tsv -t test/trait.tsv -n test/tree.nwk -o test_res.tsv; \n'+ \
    '[With non-binary gene table] python Evolink.py -g test/gene_CN.tsv -c -t test/trait.tsv -n test/tree.nwk -o test_res.tsv.')
    
    # Essential Input
    parser.add_argument('-g', '--gene', help='Tab-delimited gene presence/absence or copy number table. Columns are gene families, while rows are tip names/species/genomes in the phylogenetic tree. If copy number table is provided, please use -c option so that it will be internally converted to binary table. Presence=1, Absence=0.', required=True, dest='gene_table', metavar='FILE_PATH')
    parser.add_argument('-t', '--trait', help='Two-column (so far only one trait a time) tab-delimited trait presence/absence table. The first column is tip names and the second column is the presence/absence of this trait on the tips/species/genomes. Presence=1, Absence=0.',required=True, dest='trait_table', metavar='FILE_PATH')
    parser.add_argument('-n', '--phylogeny', help='A phylogentic tree in newick format. The tip names should be the same in the gene table and trait table.',required=True, dest='tree_file', metavar='FILE_PATH')

    # Optional Input
    parser.add_argument('-c', '--copy_number', help='The given gene table stores numbers (e.g. gene copy numbers) instead of presence/absence binary values. [Default: True]', action='store_true', required=False, default=False, dest='CN')
    # parser.add_argument('-p', '--permutation_times', help='Need to do permutation and set the permuation times [Range: 0, 10, 100, 10000, 10000]. Default is 0 and no permutation is performed.', required=False, choices=[10, 100, 1000, 10000], default=0, dest='permutation_times', metavar='INT')
    # parser.add_argument('-@', '--threads', help='Threads to use [Default: 1]', required=False, default=4, dest='threads', metavar='INT')

    # Output
    parser.add_argument('-o', '--output', help='Output file', required=True, dest='output', metavar='FILE_PATH')
    
    args = parser.parse_args()
    pipeline(args)
