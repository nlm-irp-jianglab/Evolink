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

def simulate_phen(tree_file, phen_file, perm_times, seed=1):

    rcode = '''
    suppressPackageStartupMessages(library(geiger))
    set.seed('''+str(seed)+''')
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
    # all_threads = args.threads
    threads = args.threads
    seed = args.seed
    threshold = args.threshold
    alpha = args.alpha
    output = args.output

    # ### Set theads ###
    # if all_threads <= 4:
    #     faith_pd_threads = all_threads - 1
    # elif 4 < all_threads <= 8:
    #     faith_pd_threads = 4
    # else:
    #     faith_pd_threads = 8
    # threads = all_threads - faith_pd_threads
    # os.putenv("OMP_NUM_THREADS", str(faith_pd_threads))
    # print(faith_pd_threads, threads)

    ### Read trait data ###
    print("[",datetime.datetime.now(),"]","Read trait data")
    trait_df = pd.read_csv(trait_table, sep='\t',index_col=0)
    species_list = trait_df.index.values
    trait_matrix = trait_df.values

    ### Read gene data ###
    print("[",datetime.datetime.now(),"]","Read gene data")
    gene_df = pd.read_csv(gene_table, sep='\t',index_col=0)
    gene_df = gene_df[species_list]
    if CN:
        gene_df.values = np.where(gene_df.values > 0, 1, 0)
    gene_list = gene_df.index.values
    gene_species_list = gene_df.columns.values
    gene_matrix = gene_df.values

    # test if the species in gene_table is the same to those in trait list
    if set(gene_species_list) == set(species_list):
        pass
    else:
        raise ValueError("Dimmension of genotype file do not match that of trait file")

    ### Calculate Evolink index ###
    print("[",datetime.datetime.now(),"]","Calculate Evolink index")
    df = Evolink_calculation(trait_matrix, gene_matrix, 
                             species_list, gene_list, tree_file)

    ### Perform permutation test ###
    if permutation_times > 0:
        # default labels for all genes
        df["label"] = "ns"
        df["higher_pval"] = np.nan
        df["lower_pval"] = np.nan
        df["higher_pval.adj"] = np.nan
        df["lower_pval.adj"] = np.nan

        print("[",datetime.datetime.now(),"]","Filter genes with", threshold, "% threshold")
        df['Evolink_index.abs'] = df['Evolink_index'].abs()
        # top threshold % fraction genes of abs(Evolink_index)
        kept_genes = df.nlargest(int(len(df)*threshold),'Evolink_index.abs').index
        kept_gene_df = gene_df.loc[kept_genes]
        kept_genes_matrix = kept_gene_df.values
        kept_gene_list = kept_gene_df.index.values
        print("Keep", kept_gene_df.shape[0],"/",df.shape[0], "genes")

        from multiprocessing import get_context
        from functools import partial
        import statsmodels.stats.multitest as smm

        print("[",datetime.datetime.now(),"]","Perform permutation test")

        simfun=partial( Evolink_calculation,
            gene_matrix = kept_genes_matrix,
            species_list = species_list,
            gene_list = kept_gene_list,
            tree_file = tree_file)

        perm_df = simulate_phen(tree_file, trait_table, permutation_times, seed)

        print("[",datetime.datetime.now(),"]","Simulate", permutation_times, "times")
        trait_matrix_input = []
        perm_df = perm_df.reindex(species_list)
        for (columnName, columnData) in perm_df.iteritems():
            perm_trait_mat = columnData.values
            trait_matrix_input.append(perm_trait_mat)

        # get_context("spawn") might solve the occupation of multiprocesses called by faith_pd in the unifrac module
        with get_context("spawn").Pool(processes=threads) as pool:
            # tqdm is used to show progress bar
            sim_list_arr = list(tqdm(pool.imap_unordered(simfun, trait_matrix_input), total=len(trait_matrix_input)))
            pool.close()
            pool.join()
        sim_res= np.stack(sim_list_arr,axis=0)[:,:,1].T

        evolink_value = df.loc[kept_genes, "Evolink_index"].values
        repeats_array = np.transpose([evolink_value] * permutation_times) # an array of len(evolink_value)*permutation_times length
        
        lower_pval = (np.sum(np.less(repeats_array, sim_res),axis=1)+1)/(permutation_times+1)
        higher_pval = (np.sum(np.greater(repeats_array, sim_res),axis=1)+1)/(permutation_times+1)
        
        df.loc[kept_genes, "higher_pval"] = higher_pval
        df.loc[kept_genes, "lower_pval"] = lower_pval
        df.loc[kept_genes, "higher_pval.adj"] = smm.multipletests(df.loc[kept_genes, "higher_pval"].values, method="fdr_bh")[1]
        df.loc[kept_genes, "lower_pval.adj"] = smm.multipletests(df.loc[kept_genes, "lower_pval"].values, method="fdr_bh")[1]
        df.loc[((df["lower_pval"] < 0.05) & (df["lower_pval.adj"] < alpha)) | (df["higher_pval"] < 0.05) & (df["higher_pval.adj"] < alpha), "label"] = "sig"

    ### Output ###
    print("[",datetime.datetime.now(),"]","Output result")
    df.to_csv(output, sep="\t", na_rep='NaN')

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
    parser.add_argument('-p', '--permutation_times', help='Need to do permutation and set the permuation times [Range: 0, 10, 100, 10000, 10000]. Default is 0 and no permutation is performed.', required=False, type=int, default=0, dest='permutation_times', metavar='PERM_TIMES') #choices=[0, 10, 100, 1000, 10000]
    parser.add_argument('-@', '--threads', help='Threads for permutation test [Default: 4]', required=False, default=4, type=int, dest='threads', metavar='THREADS')
    parser.add_argument('-s', '--seed', help='Set seed for simulation for reproducibility of the results [Default: 1]', required=False, default=1, type=int, dest='seed', metavar='SEED')
    parser.add_argument('-f', '--threshold', help='Evolink index top percentage cutoff to filter in genes for permutation tests [Range: 0-1; Default: 0.02]', required=False, default=0.02, type=float, dest='threshold', metavar='THRESHOLD')
    parser.add_argument('-a', '--alpha', help='Adjusted p-value cutoff [Range: 0-1; Default: 0.05]', required=False, default=0.05, type=float, dest='alpha', metavar='ALPHA')
    
    # Output
    parser.add_argument('-o', '--output', help='Output file', required=True, dest='output', metavar='OUTPUT')
    
    args = parser.parse_args()
    pipeline(args)
