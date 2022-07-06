#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os, shutil, sys, tempfile, datetime, subprocess
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from unifrac import faith_pd
from biom.util import biom_open
from biom.table import Table
from tqdm import tqdm

script_dir = os.path.abspath(os.path.dirname( __file__ ))
pvalue2zscore = {0.80:1.282, 0.85:1.44, 0.9:1.645, 0.95:1.960, 0.99:2.576, 0.995:2.807, 0.997:3.0, 0.999:3.291}

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
    tree_bi=multi2di(tree) # force tree to be bifurcate
    tree_bi$edge.length[tree_bi$edge.length==0]=min(tree$edge.length) # not allow br length is zero
    sim_data = simTrait(tree_bi, phen, perm_times='''+str(perm_times)+''')
    '''
    sim_data = robjects.r(rcode)

    with localconverter(robjects.default_converter + pandas2ri.converter):
        sim_df = robjects.conversion.rpy2py(sim_data)

    return(sim_df)

def sig_genes(df, alpha=0.01, p_alpha=0.01, permutation=False):
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index':'orthoID'})
    if permutation:
        df['significance'] = np.where(((df["pvalue"] < alpha) & (df["permutation_pvalue"] < p_alpha)), "sig", pd.NA)
    else:
        df['significance'] = np.where((df["pvalue"] < alpha), "sig", pd.NA)
    return(df)

def extract_list(lst, index_list):
    return(",".join([str(lst[i]) for i in [0]+index_list]))

def iTOL_data(trait_df, gene_df, pos_genes, neg_genes):
    # set field annotation
    field_shapes = ",".join(['2'] + ['1']*len(pos_genes) + ['1']*len(neg_genes))
    field_labels = ",".join(["Phenotype"] + pos_genes + neg_genes)
    pos_palette = ['#93003a', '#a62045', '#b73651', '#c84b5e', '#d7606d', '#e4757d', '#ef8a8f', '#f8a1a2', '#ffb7b7', '#ffd0cf']
    len_pos_pal = len(pos_palette)
    neg_palette = ['#00429d', '#2255a0', '#3866a5', '#4d77ac', '#6188b4', '#7699bd', '#8baac6', '#a1bbd1', '#b7ccdc', '#cfdde7']
    len_neg_pal = len(neg_palette)
    field_colors = ",".join(['#FFB844'] + pos_palette*int(len(pos_genes)/len_pos_pal) + pos_palette[:len(pos_genes)%len_pos_pal] \
                   + neg_palette*int(len(neg_genes)/len_neg_pal) + neg_palette[:len(neg_genes)%len_neg_pal]) \
    # set legend annotation
    legend_title = "Annotated Tree"
    legend_position_x = '100'
    legend_position_y = '100'
    legend_shapes = field_shapes
    legend_labels = field_labels
    legend_colors = field_colors
    subset_gene_df = gene_df.loc[pos_genes + neg_genes]

    comb_df = pd.concat([trait_df,subset_gene_df.T], axis=1)
    data = list(comb_df.to_records(index=True))
    data = [tuple(e) for e in data]
    return(data, legend_title, legend_position_x, legend_position_y, 
           legend_shapes, legend_labels, legend_colors,
           field_shapes, field_labels, field_colors)

def iTOL_input(tree_file, trait_df, gene_df, df, prefix, pos_top=5, neg_top=5):
    # display_mode:1=rectangular, 2=circular.
    import zipfile
    
    itol_tree_path = os.path.join(prefix, "input.tree")
    shutil.copy(tree_file, itol_tree_path)
    zip_path = os.path.join(prefix, "Evolink_itol_input.zip")
    with zipfile.ZipFile(zip_path, mode="w") as archive:
        archive.write(itol_tree_path, os.path.basename(itol_tree_path))

    pos_genes = df[(df["significance"]=="sig") & (df["Evolink_index"] > 0)].sort_values(by='Evolink_index', ascending=False).head(pos_top)["orthoID"].to_list()
    neg_genes = df[(df["significance"]=="sig") & (df["Evolink_index"] < 0)].sort_values(by='Evolink_index').head(neg_top)["orthoID"].to_list()

    sys.path.append(script_dir)
    from iTOL import TOL
    t = TOL(zfile=zip_path, wd=prefix)
    # for each tip, it will have 1 phenotype + pos_top genotypes + neg_top genotypes
    data, legend_title, legend_position_x, legend_position_y, legend_shapes, legend_labels, legend_colors, field_shapes, field_labels, field_colors = iTOL_data(trait_df, gene_df, pos_genes, neg_genes)
    t.binary(data, separator='comma', dataset_label='binary', 
             field_shapes=field_shapes, field_labels=field_labels,
             field_colors=field_colors, legend_title=legend_title, legend_position_x=legend_position_x, 
             legend_position_y=legend_position_y, legend_shapes=legend_shapes, 
             legend_labels=legend_labels, legend_colors=legend_colors)

    annot_path = os.path.join(prefix, "binary.txt")
    with zipfile.ZipFile(zip_path, mode="a") as archive:
        archive.write(annot_path, os.path.basename(annot_path))
    
    print("[",datetime.datetime.now(),"]", "iTOL inputs generated at", zip_path, flush=True)

def plot_fig(gene_table, trait_table, tree_file, outdir, cutoff, pos_top=5, neg_top=5, display_mode=1):
    outfile = os.path.join(outdir, "result.tsv")
    script_dir = os.path.abspath(os.path.dirname( __file__ ))

    cmd = "Rscript --vanilla {0}/Evolink_plot.R -g {1} -t {2} -n {3} -r {4} -c {5} -a {6} -b {7} -m {8} -o {9}".format(
        script_dir,
        gene_table, trait_table, tree_file, outfile, 
        cutoff, pos_top, neg_top, display_mode, outdir
    )
    print("Plot command line:", cmd, flush=True)
    subprocess.call(cmd, shell=True)

def pipeline(args):

    gene_table = args.gene_table
    trait_table = args.trait_table
    tree_file = args.tree_file
    CN = args.CN
    alpha = args.alpha
    p_threshold = args.p_threshold
    e_threshold = args.e_threshold
    permutation_times = args.permutation_times
    threads = args.threads
    seed = args.seed
    p_alpha = args.p_alpha
    # multitest_correction = args.multitest_correction
    plot = args.plot
    top_genes = args.top_genes
    display_mode = args.display_mode
    force = args.force
    output = args.output

    if os.path.exists(output) and os.path.isdir(output):
        if force:
            shutil.rmtree(output)
            os.makedirs(output)
        else:
            sys.exit("Error: {} already exists. To force overwrite it use -f or --force.".format(output))
    else:
        os.makedirs(output)

    pos_top, neg_top = [int(i) for i in top_genes.split(",")]

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
    
    # default pvalues for all genes
    df["pvalue"] = pd.NA
    # df["adjusted_pvalue"] = np.nan

    pos_ct=np.sum(trait_matrix, axis=0)[0]
    all_ct = trait_matrix.shape[0]
    minority_ct = min(pos_ct, all_ct-pos_ct)
    # genes only existing in 25% of minority will still be kept
    rare_cutoff = 2*(minority_ct*0.25/all_ct)-1 
    majority_ct = all_ct - minority_ct
    prevalence_cutoff = 2*(majority_ct/all_ct)-1

    ### Filter genes ###
    print("[",datetime.datetime.now(),"]","Filter genes", flush=True)
    kept_genes = df.loc[((df['Prevalence_index'] >= np.min([-1*p_threshold, rare_cutoff])) & (df['Prevalence_index'] <= np.max([p_threshold, prevalence_cutoff]))),].index
    kept_gene_df = gene_df.loc[kept_genes]
    kept_gene_matrix = kept_gene_df.values
    kept_gene_list = kept_gene_df.index.values
    print("[",datetime.datetime.now(),"]","Keep", str(kept_gene_df.shape[0])+"/"+str(df.shape[0]), "genes", flush=True)

    ### Get basic p values ####
    from statsmodels.distributions.empirical_distribution import ECDF
    ecdf = ECDF(np.abs(df.loc[kept_genes, "Evolink_index"]))
    pval = 1 - ecdf(np.abs(df.loc[kept_genes, "Evolink_index"]))
    df.loc[kept_genes, "pvalue"] = pval
    print("[",datetime.datetime.now(),"]","Get Evolink p values", flush=True)

    ### Get Evolink index threshold ###
    if e_threshold:
        thresh = e_threshold
        print("[",datetime.datetime.now(),"]","Get Evolink index threshold from user:", thresh, flush=True)
        # alpha will be updated based on threshold provided by the user
        alpha = 1 - ecdf(thresh)
        print("[",datetime.datetime.now(),"]", "According pvalue threshold (alpha) is", alpha, flush=True)
    else:
        thresh = np.quantile(np.abs(df.loc[kept_genes, "Evolink_index"]), q=1-alpha)
        print("[",datetime.datetime.now(),"]","Get Evolink index threshold:", thresh, flush=True)

    if permutation_times > 0:
        permutation = True
        print("[",datetime.datetime.now(),"]","Perform permutation test", flush=True)
        
        from multiprocessing import get_context
        from functools import partial
        # import statsmodels.stats.multitest as smm

        print("[",datetime.datetime.now(),"]","Simulate", permutation_times, "times", flush=True)
        simfun = partial( Evolink_calculation,
            gene_matrix = kept_gene_matrix,
            species_list = species_list,
            gene_list = kept_gene_list,
            tree_file = tree_file)

        perm_df = simulate_phen(tree_file, trait_table, permutation_times, seed)

        trait_matrix_input = []
        perm_df = perm_df.reindex(species_list)
        for (_, columnData) in perm_df.iteritems():
            perm_trait_mat = columnData.values
            trait_matrix_input.append(perm_trait_mat)

        # get_context("spawn") might solve the occupation of multiprocesses called by faith_pd in the unifrac module
        with get_context("spawn").Pool(processes=threads) as pool:
            # tqdm is used to show progress bar
            sim_list_arr = list(tqdm(pool.imap_unordered(simfun, trait_matrix_input), total=len(trait_matrix_input)))
            pool.close()
            pool.join()
        sim_res= np.stack(sim_list_arr,axis=0)[:,:,1].T #1 for the 2nd column of result table
        print()

        # sim_thresh = np.quantile(np.abs(sim_res), q=1-alpha) # thresh for Evolink index using simulation data
        
        evolink_value = df.loc[kept_genes, "Evolink_index"].values
        repeats_array = np.transpose([evolink_value] * permutation_times) # an array of len(evolink_value)*permutation_times length

        pval = (np.sum(np.less(abs(repeats_array), abs(sim_res)),axis=1)+1)/(permutation_times+1)
        df.loc[kept_genes, "permutation_pvalue"] = pval
        # df.loc[kept_genes, "adjusted_pvalue"] = smm.multipletests(df.loc[kept_genes, "pvalue"].values, method=multitest_correction)[1]
        print("[",datetime.datetime.now(),"]","Get permutation p values", flush=True)

    else:
        permutation = False

    res_df = sig_genes(df, alpha, p_alpha, permutation)
    sig_gene_ct = res_df[res_df["significance"]=="sig"].shape[0]
    print("[",datetime.datetime.now(),"]","Find", sig_gene_ct, "significant genes", flush=True)
    
    ### Output ###
    print("[",datetime.datetime.now(),"]","Output result", flush=True)
    outfile = os.path.join(output, "result.tsv")
    res_df.to_csv(outfile, sep="\t", index=False, na_rep='NA')

    ### Plot ####
    if plot:
        print("[",datetime.datetime.now(),"]","Generate figures", flush=True)
        iTOL_input(tree_file, trait_df, gene_df, res_df, output, pos_top, neg_top)
        plot_fig(gene_table, trait_table, tree_file, output, alpha, pos_top, neg_top, display_mode)

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Evolink is designed to find gene families associated with trait by explicitly using phylogeny information.')
    
    # Essential Input
    parser.add_argument('-g', '--genotype', help='Tab-delimited gene presence/absence or copy number table. Columns are gene families, while rows are tip names/species/genomes in the phylogenetic tree. If copy number table is provided, please use -c option so that it will be internally converted to binary table. Presence=1, Absence=0.', required=True, dest='gene_table', metavar='GENE_TABLE')
    parser.add_argument('-t', '--phenotype', help='Two-column (so far only one trait is allowed each time) tab-delimited trait presence/absence table. The first column is tip names and the second column is the presence/absence of this trait on the tips/species/genomes. Presence=1, Absence=0.',required=True, dest='trait_table', metavar='TRAIT_TABLE')
    parser.add_argument('-n', '--phylogeny', help='A phylogentic tree in newick format. The tip names should be the same in the gene table and trait table.',required=True, dest='tree_file', metavar='TREE')

    # Optional Input
    parser.add_argument('-c', '--copy_number', help='The given gene table stores numbers (e.g. gene copy numbers) instead of presence/absence binary values. [Default: True]', action='store_true', required=False, default=False, dest='CN')
    parser.add_argument('-a', '--alpha', help='Pvalue threshold [Default:0.01]', required=False, type=float, default=0.01, dest='alpha', metavar='ALPHA')
    parser.add_argument('-p', '--p_threshold', help='Absolute Prevalence index threshold to filter genes and get Evolink index distribution [Range: 0-1; Default: 0.9]', required=False, default=0.9, type=float, dest='p_threshold', metavar='THRESHOLD')
    parser.add_argument('-e', '--e_threshold', help='Absolute Evolink index threshold to select significant genes. Notice: P-value cutoff (alpha) will be updated based on this option [Range: 0-1; Default: NULL]', required=False, default=None, type=float, dest='e_threshold', metavar='THRESHOLD')
    
    # Simulation Options
    parser.add_argument('-s', '--simulation_times', help='Need to permutation test and set the simulation times [Range: 0-10000]. Default is 0 and no permutation is performed.', required=False, type=int, default=0, dest='permutation_times', metavar='PERM_TIMES')
    parser.add_argument('-@', '--threads', help='Threads for permutation test [Default: 4]', required=False, default=4, type=int, dest='threads', metavar='THREADS')
    parser.add_argument('-r', '--seed', help='Set seed for simulation for reproducibility of the results [Default: 1]', required=False, default=1, type=int, dest='seed', metavar='SEED')
    parser.add_argument('-b', '--permutation_alpha', help='Permutation pvalue threshold [Default:0.01]', required=False, type=float, default=0.01, dest='p_alpha', metavar='PERMUTATION_ALPHA')
    # parser.add_argument('-m', '--multitest_correction', help='Multitest correction [Choices: bonferroni, fdr_bh, holm, hommel; Default: bonferroni]', choices = ["bonferroni", "fdr_bh", "holm", "hommel"], required=False, type=str, default="bonferroni", dest='multitest_correction', metavar='MT_CORR')
    
    # Plot Options
    parser.add_argument('-v', '--visualization', help='Whether to generate plots', action='store_true', required=False, default=False, dest='plot')
    parser.add_argument('-N', '--top_genes', help='Top positively and negatively associated genes mapped to tree. [Default: 5,5 for top 5 pos genes and top 5 neg genes.]', required=False, default='5,5', type=str, dest='top_genes')
    parser.add_argument('-d', '--display-mode', help='Tree display mode. [1: circular, 2: rectangular; Default: 1]', type=int, choices=[1, 2], required=False, default=1, dest='display_mode')
    parser.add_argument('-f', '--force', help='Force to overwrite output folder. [Default: False]', action='store_true', required=False, default=False, dest='force')
    
    # Output
    parser.add_argument('-o', '--output', help='output directory', required=True, dest='output', metavar='OUTPUT')
    
    args = parser.parse_args()
    pipeline(args)

