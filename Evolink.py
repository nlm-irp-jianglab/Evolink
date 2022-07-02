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
# from tqdm import tqdm

script_dir = os.path.abspath(os.path.dirname( __file__ ))

alpha2CI = {0.80:1.282, 0.85:1.44, 0.9:1.645, 0.95:1.960, 0.99:2.576, 0.995:2.807, 0.997:3.0, 0.999:3.291}


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

    result = np.vstack((G1T1_T1, G1T0_T0, G0T0_T0, G0T1_T1, Prevalence_index, Evolink_index)).T
    df = pd.DataFrame(result, columns=["G1T1_T1", "G1T0_T0", "G0T0_T0", "G0T1_T1", "Prevalence_index", "Evolink_index"], index=gene_list)
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

def sig_genes(df, prevalence_cutoff, rare_cutoff, p_threshold, e_threshold, alpha):
    with localconverter(robjects.default_converter + pandas2ri.converter):
        r_df = robjects.conversion.py2rpy(df)
    robjects.globalenv['r_df'] = r_df
    robjects.globalenv['prevalence_cutoff'] = float(prevalence_cutoff)
    robjects.globalenv['rare_cutoff'] = float(rare_cutoff)
    robjects.globalenv['p_threshold'] = float(p_threshold)
    robjects.globalenv['e_threshold'] = float(e_threshold)
    robjects.globalenv['CI'] = alpha2CI[alpha]
    rcode = '''
    suppressPackageStartupMessages({
        suppressWarnings(library(tidyverse))
    })
    r_df = r_df %>% rownames_to_column(var = "orthoID")
    dt = r_df %>% filter(Prevalence_index >= min(-1*p_threshold, rare_cutoff), Prevalence_index <= max(p_threshold, prevalence_cutoff))
    std = sd(dt$Evolink_index)
    dt = dt %>% mutate(z_score=(Evolink_index-0)/std, significance=ifelse(abs(z_score)>=CI & abs(Evolink_index) >= e_threshold, "sig", NA))
    r_df = r_df %>% select(orthoID, Prevalence_index, Evolink_index) %>% 
                left_join(dt, by = c("orthoID", "Prevalence_index", "Evolink_index")) %>% 
                mutate(z_score=(Evolink_index-0)/std) %>% 
                select(orthoID, Prevalence_index, Evolink_index, significance, z_score)
    '''
    r_df = robjects.r(rcode)
    with localconverter(robjects.default_converter + pandas2ri.converter):
        res_df = robjects.conversion.rpy2py(r_df)
    return(res_df)

def extract_list(lst, index_list):
    return(",".join([str(lst[i]) for i in [0]+index_list]))

def iTOL_data(trait_df, gene_df, pos_genes, neg_genes):
    # set field annotation
    field_shapes = ",".join(['2'] + ['1']*len(pos_genes) + ['1']*len(neg_genes))
    field_labels = ",".join(["Trait"] + pos_genes + neg_genes)
    pos_palette = ['#93003a', '#a62045', '#b73651', '#c84b5e', '#d7606d', '#e4757d', '#ef8a8f', '#f8a1a2', '#ffb7b7', '#ffd0cf']
    len_pos_pal = len(pos_palette)
    neg_palette = ['#00429d', '#2255a0', '#3866a5', '#4d77ac', '#6188b4', '#7699bd', '#8baac6', '#a1bbd1', '#b7ccdc', '#cfdde7']
    len_neg_pal = len(neg_palette)
    field_colors = ",".join(['#FFB844'] + pos_palette*int(len(pos_genes)/len_pos_pal) + pos_palette[:len(pos_genes)%len_pos_pal] \
                   + neg_palette*int(len(neg_genes)/len_neg_pal) + neg_palette[:len(neg_genes)%len_neg_pal]) \
    # set legend annotation
    legend_title = "Annotated Tree"
    legend_position_x = '0'
    legend_position_y = '1000'
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

def plot_fig(gene_table, trait_table, tree_file, outdir, CI=3, pos_top=5, neg_top=5, display_mode=1):
    outfile = os.path.join(outdir, "result.tsv")
    script_dir = os.path.abspath(os.path.dirname( __file__ ))

    cmd = "Rscript --vanilla {0}/Evolink_plot.R -g {1} -t {2} -n {3} -r {4} -c {5} -a {6} -b {7} -m {8} -o {9}".format(
        script_dir,
        gene_table, trait_table, tree_file, outfile, 
        CI, pos_top, neg_top, display_mode, outdir
    )
    print("Plot command line:", cmd, flush=True)
    subprocess.call(cmd, shell=True)

def pipeline(args):

    gene_table = args.gene_table
    trait_table = args.trait_table
    tree_file = args.tree_file
    CN = args.CN
    p_threshold = args.p_threshold
    e_threshold = args.e_threshold
    alpha = args.alpha
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

    ### Print command line ###
    cn = "-c " if CN else ""
    mode = "-m " if display_mode else ""
    vis = "-v " if plot else ""
    f = "-f " if force else ""
    cmd = "python Evolink.py -g {0} -t {1} -n {2} {3}-p {4} -e {5} -a {6} {7}-N {8} {9}-o {10} {11}".format(
        gene_table, trait_table,
        tree_file, cn,
        p_threshold, e_threshold,
        alpha, vis, top_genes, 
        mode, output, f
    )
    print("Main command line:", cmd, flush=True)

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
    pos_ct=np.sum(trait_matrix, axis=0)[0]
    all_ct = trait_matrix.shape[0]
    minority_ct = min(pos_ct, all_ct-pos_ct)
    # genes only existing in 25% of minority will still be kept
    rare_cutoff = 2*(minority_ct*0.25/all_ct)-1 
    majority_ct = all_ct - minority_ct
    prevalence_cutoff = 2*(majority_ct/all_ct)-1
    res_df = sig_genes(df, prevalence_cutoff, rare_cutoff, p_threshold, e_threshold, alpha)
    CI = alpha2CI[alpha]

    ### Output ###
    print("[",datetime.datetime.now(),"]","Output result", flush=True)
    outfile = os.path.join(output, "result.tsv")
    res_df.to_csv(outfile, sep="\t", index=False, na_rep='NA')

    ### Plot ####
    if plot:
        print("[",datetime.datetime.now(),"]","Generate figures", flush=True)
        iTOL_input(tree_file, trait_df, gene_df, res_df, output, pos_top, neg_top)
        plot_fig(gene_table, trait_table, tree_file, output, CI, pos_top, neg_top, display_mode)

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Evolink is designed to find gene families associated with given trait with the help of phylogeny information.')
    
    # Essential Input
    parser.add_argument('-g', '--genotype', help='Tab-delimited gene presence/absence or copy number table. Columns are gene families, while rows are tip names/species/genomes in the phylogenetic tree. If copy number table is provided, please use -c option so that it will be internally converted to binary table. Presence=1, Absence=0.', required=True, dest='gene_table', metavar='GENE_TABLE')
    parser.add_argument('-t', '--phenotype', help='Two-column (so far only one trait is allowed each time) tab-delimited trait presence/absence table. The first column is tip names and the second column is the presence/absence of this trait on the tips/species/genomes. Presence=1, Absence=0.',required=True, dest='trait_table', metavar='TRAIT_TABLE')
    parser.add_argument('-n', '--phylogeny', help='A phylogentic tree in newick format. The tip names should be the same in the gene table and trait table.',required=True, dest='tree_file', metavar='TREE')

    # Optional Input
    parser.add_argument('-c', '--copy_number', help='The given gene table stores numbers (e.g. gene copy numbers) instead of presence/absence binary values. [Default: True]', action='store_true', required=False, default=False, dest='CN')
    parser.add_argument('-p', '--p_threshold', help='Absolute Prevalence index threshold to filter genes and get background distribution [Range: 0-1; Default: 0.9]', required=False, default=0.9, type=float, dest='p_threshold', metavar='THRESHOLD')
    parser.add_argument('-e', '--e_threshold', help='Absolute Evolink index threshold to select significant genes [Range: 0-1; Default: 0.1]', required=False, default=0.1, type=float, dest='e_threshold', metavar='THRESHOLD')
    parser.add_argument('-a', '--alpha', help='Tail of area cutoff [Range: 0.8, 0.85, 0.90, 0.95, 0.99, 0.995, 0.997, 0.999; Default: 0.997]', choices=[0.8, 0.85, 0.90, 0.95, 0.99, 0.995, 0.997, 0.999], required=False, default=0.997, type=float, dest='alpha', metavar='ALPHA')
    parser.add_argument('-v', '--visualization', help='Whether to generate plots', action='store_true', required=False, default=False, dest='plot')
    parser.add_argument('-N', '--top_genes', help='Top positively and negatively associated genes mapped to tree. [Default: 5,5 for top 5 pos genes and top 5 neg genes.]', required=False, default='5,5', type=str, dest='top_genes')
    parser.add_argument('-m', '--display-mode', help='Tree display mode. [1: circular, 2: rectangular; Default: 1]', type=int, choices=[1, 2], required=False, default=1, dest='display_mode')
    parser.add_argument('-f', '--force', help='Force to overwrite output folder. [Default: False]', action='store_true', required=False, default=False, dest='force')
    
    # Output
    parser.add_argument('-o', '--output', help='output directory', required=True, dest='output', metavar='OUTPUT')
    
    args = parser.parse_args()
    pipeline(args)
