#!/usr/bin/env python
# -*- coding: utf-8 -*-

from tqdm import tqdm
import pandas as pd
import numpy as np
import os, tempfile

def Evolink_calculation(trait_matrix, species_list, gene_species_list, gene_list, gene_matrix, tree_file):

    T1_index = np.where(trait_matrix == 1)[0]
    T0_index = np.where(trait_matrix == 0)[0]
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

    result = np.vstack((Prevalence_index, Evolink_index)).T
    df = pd.DataFrame(result, columns=["Prevalence_index", "Evolink_index"], index=gene_list)

    return df

def matrix_faith_pd( matrix, species_name, gene_list,tree_path):
    from unifrac import faith_pd
    fd, path = tempfile.mkstemp(suffix=".biom")
    matrix2hdf5(matrix, species_name, gene_list, path)
    obs = faith_pd(path,tree_path)
    os.close(fd)
    os.unlink(path)
    return obs

def compute_four_matrix(gene_species_list, species_list, T1_species,T0_species, gene_matrix, trait_matrix):
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

def create_permuation_matrix(list1, list2):
    matrix= np.zeros((len(list1), len(list2)))
    for element in list2:
        idx1 = list1.index(element)
        idx2 = list2.index(element)
        matrix[idx1,idx2] = 1
    return matrix

def simulate_phen(tree_file, phen_file, perm_times):
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter

    rcode = '''
    suppressPackageStartupMessages(library(geiger))
    cn2bi<-function(x, n){
        cutoff=sort(x, decreasing = TRUE)[n]
        return(as.integer(x >= cutoff))
    }
    simTrait=function(tree, phen, perm_times=10000){
        phen_pos_n=sum(phen)
        rateMat=ratematrix(tree, phen)
        sim_res=sim.char(tree, rateMat, nsim = perm_times)
        sim_table=data.frame(sim_res)
        colnames(sim_table)=paste0("N", 1:perm_times)
        bi_sim_table=data.frame(apply(sim_table, 2, cn2bi, n=phen_pos_n))
        rownames(bi_sim_table)=rownames(sim_res)
        return(bi_sim_table)
    }
    tree=read.tree("'''+tree_file+'''")
    phen=read.table("'''+phen_file+'''", header=T, row.names=1)
    sim_data = simTrait(tree, phen, perm_times='''+str(perm_times)+''')
    '''

    sim_data = robjects.r(rcode)

    with localconverter(robjects.default_converter + pandas2ri.converter):
        sim_df = robjects.conversion.rpy2py(sim_data)

    return(sim_df)

def pipeline(args):

    gene_table = args.gene_table
    trait_table = args.trait_table
    tree_file = args.tree_file
    CN = args.CN
    permutation_times = args.permutation_times
    threads = args.threads
    output = args.output

    ### Read trait data ###
    print("Read trait data")
    trait_df = pd.read_csv(trait_table, sep='\t',index_col=0)
    species_list = trait_df.index.values.tolist()
    trait_matrix = np.ravel(trait_df.values)

    ### Read gene data ###
    print("Read gene data")
    gene_df = pd.read_csv(gene_table, sep='\t',index_col=0)
    if CN:
        gene_df.values = np.where(gene_df.values > 0, 1, 0)
    gene_list = gene_df.index.values
    gene_species_list = gene_df.columns.values.tolist()
    gene_matrix = gene_df.values

    # test if the species in gene_table is the same to those in trait list
    if set(gene_species_list) == set(species_list):
        pass
    else:
        raise ValueError("dimmension of genotype file do not match that of trait file")

    ### Perform permutation test ###
    perm_done = 0
    if permutation_times > 0:
        print("Perform permutation test")
        from multiprocessing import Pool
        from functools import partial

        simfun=partial( Evolink_calculation,
            species_list = species_list,
            gene_species_list=gene_species_list,
            gene_list=gene_list,
            gene_matrix=gene_matrix,
            tree_file=tree_file)

        perm_df = simulate_phen(tree_file, trait_table, permutation_times)

        trait_matrix_input = []
        perm_df = perm_df.reindex(species_list)
        for (columnName, columnData) in perm_df.iteritems():
            perm_trait_mat = columnData.values
            trait_matrix_input.append(perm_trait_mat)

        print("Multiprocess simulate", permutation_times, "times")
        with Pool(processes=threads) as pool:
            # tqdm is used to show progress bar
            sim_list_arr = list(tqdm(pool.imap_unordered(simfun, trait_matrix_input), total=len(trait_matrix_input)))
            sim_res= np.stack(sim_list_arr,axis=0)[:,:,1].T
        pool.close()
        pool.join()
        perm_done = 1

    print("Calculate Evolink index")
    df = Evolink_calculation(trait_matrix, species_list, gene_species_list, gene_list,
                             gene_matrix, tree_file)

    evolink_value = df["Evolink_index"].values
    if perm_done:
        repeats_array = np.transpose([evolink_value] * permutation_times) # an array of len(evolink_value)*permutation_times length
        lower_pval = (np.sum(np.less(repeats_array, sim_res),axis=1)+1)/(permutation_times+1)
        higher_pval = (np.sum(np.greater(repeats_array, sim_res),axis=1)+1)/(permutation_times+1)
        df["higher_pval"] = higher_pval
        df["lower_pval"] = lower_pval

    print("Output result")
    df.to_csv(output, sep="\t", na_rep='NaN')

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
    parser.add_argument('-p', '--permutation_times', help='Need to do permutation and set the permuation times [Range: 0, 10, 100, 10000, 10000]. Default is 0 and no permutation is performed.', required=False, type=int, default=0, dest='permutation_times', metavar='INT') #choices=[0, 10, 100, 1000, 10000]
    parser.add_argument('-@', '--threads', help='Threads to use [Default: 1]', required=False, default=1, type=int, dest='threads', metavar='INT')

    # Output
    parser.add_argument('-o', '--output', help='Output file', required=True, dest='output', metavar='FILE_PATH')
    
    args = parser.parse_args()
    pipeline(args)
