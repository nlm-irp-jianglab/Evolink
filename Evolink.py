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

def Evolink_calculation(gene_matrix, trait_matrix, species_list, gene_list, tree_file):

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
    t = read.tree("'''+tree_file+'''")
    sum(t$edge.length)
    ''')
    br_len_sum = res[0]
    return br_len_sum

def read_tree(tree_file):
    r_tree = robjects.r('''
    library(ape)
    tree = read.tree("'''+tree_file+'''")
    ''')
    robjects.globalenv['tree'] = r_tree

def generate_tree(tree_file, species_sub):
    fd, path = tempfile.mkstemp(suffix=".nwk")
    species = robjects.vectors.StrVector(species_sub)
    robjects.globalenv['species'] = species
    read_tree(tree_file)
    robjects.r('''
    pruned.tree<-drop.tip(tree, setdiff(tree$tip.label, species));
    write.tree(pruned.tree, "'''+path+'''")
    ''')
    br_len = sum_br_len(path)
    os.close(fd)
    return (path, br_len)

def simulate_gene(filtered_gene, sim_fold=10, seed=1, tree_file=None):
    # tree is already in robject.globalenv
    if tree_file!=None:
        read_tree(tree_file)

    with localconverter(robjects.default_converter + pandas2ri.converter):
        r_filtered_gene = robjects.conversion.py2rpy(filtered_gene)
    robjects.globalenv['fgene'] = r_filtered_gene

    rcode = '''
    source("'''+script_dir+'''/simulate_gene.R")
    sim_res = gene.sim(tree=tree, gene=fgene, sim.fold='''+str(sim_fold)+''', seed='''+str(seed)+''')
    '''

    sim_res = robjects.r(rcode) # row: genes; col: species/tips
    with localconverter(robjects.default_converter + pandas2ri.converter):
        sim_df = robjects.conversion.rpy2py(sim_res)
    return(sim_df)

def gesd_test(df, kept_genes):
    r_kept_genes = robjects.vectors.StrVector(kept_genes)
    robjects.globalenv['kept_genes'] = r_kept_genes

    with localconverter(robjects.default_converter + pandas2ri.converter):
        r_df = robjects.conversion.py2rpy(df)
    robjects.globalenv['df'] = r_df
    
    rcode = '''
    source("'''+script_dir+'''/gesd_test.R")
    x = df[which(rownames(df) %in% kept_genes), "Evolink_index"]
    suppressWarnings({out = gesdTest(x)})
    df[which(rownames(df) %in% kept_genes), "pvalue"] = out$p.value
    r_df = df
    '''
    r_df = robjects.r(rcode)
    with localconverter(robjects.default_converter + pandas2ri.converter):
        df = robjects.conversion.rpy2py(r_df)
    return(df)

### Plot functions ###
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

def plot_fig(mode, gene_table, trait_table, tree_file, outdir, cutoff=0.05, pos_top=5, neg_top=5, display_mode=1):
    # add mode
    # redefine cutoff
    outfile = os.path.join(outdir, "result.tsv")

    cmd = "Rscript --vanilla {0}/Evolink_plot.R -m {1} -g {2} -t {3} -n {4} -r {5} -c {6} -a {7} -b {8} -d {9} -o {10}".format(
        script_dir, mode, 
        gene_table, trait_table, tree_file, outfile, 
        cutoff, pos_top, neg_top, display_mode, outdir
    )
    print("Run plot command line:", cmd, flush=True)
    subprocess.call(cmd, shell=True)

def pipeline(args):

    gene_table = args.gene_table
    trait_table = args.trait_table
    tree_file = args.tree_file
    CN = args.CN
    p_threshold = args.p_threshold
    mode = args.mode
    seed = args.seed

    # Cutoff mode
    e_threshold = args.e_threshold
    
    # Resd test mode
    gesd_mc_method = args.gesd_mc_method
    gesd_pval_threshold = args.gesd_pval_threshold
    gesd_padj_threshold = args.gesd_padj_threshold

    # Permutation mode
    fold_times = args.fold_times
    perm_mc_method = args.perm_mc_method
    perm_padj_threshold = args.perm_padj_threshold
    # threads = args.threads

    # Z-score mode
    alpha = args.alpha

    # IsolationForest mode #
    outlier_score_threshold = args.outlier_score_threshold
    # slope_cutoff = args.slope_cutoff
    # threshold_interval = args.threshold_interval
    n_estimators = args.n_estimators
    max_samples = args.max_samples

    # Plot options
    plot = args.plot
    top_genes = args.top_genes
    display_mode = args.display_mode
    force = args.force
    output = args.output

    pval_threshold = 0.05

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
    df = Evolink_calculation(gene_matrix, trait_matrix,
                             species_list, gene_list, tree_file)
    
    print("[",datetime.datetime.now(),"]","Running", mode, "mode to detect significant genes", flush=True)

    if mode in ["z_score", "permutation", "gesd_test"]:
    #### For Z-score and permutation modes, do filtering first ####
        ### Filter genes ###
        print("[",datetime.datetime.now(),"]","Filter genes", flush=True)
        pos_ct=np.sum(trait_matrix, axis=0)[0]
        all_ct = trait_matrix.shape[0]
        minority_ct = min(pos_ct, all_ct-pos_ct)
        # genes only existing in 25% of minority will still be kept
        # rare_cutoff = 2*(minority_ct*0.25/all_ct)-1
        # majority_ct = all_ct - minority_ct
        # prevalence_cutoff = 2*(majority_ct/all_ct)-1

        kept_genes = df.loc[((df['Prevalence_index'] >= -1*p_threshold) & (df['Prevalence_index'] <= p_threshold)),].index
        kept_gene_df = gene_df.loc[kept_genes]
        print("[",datetime.datetime.now(),"]","Keep", str(kept_gene_df.shape[0])+"/"+str(df.shape[0]), "genes", flush=True)
        
        import statsmodels.stats.multitest as smm
        if mode == "gesd_test":
        ##### 1. Gesd Test mode end #####
            df = gesd_test(df, kept_genes)
            df["pvalue.adj"] = pd.NA
            df["significance"] = pd.NA
            
            if gesd_mc_method == "none":
                df.loc[(df["pvalue"] < gesd_pval_threshold), "significance"] = "sig"
                df["pvalue.adj"] = df["pvalue"]
                pval_threshold = gesd_pval_threshold
            else:
                df.loc[df.index.isin(kept_genes) & (df["pvalue"].notna()), "pvalue.adj"] = smm.multipletests(df.loc[df.index.isin(kept_genes) & (df["pvalue"].notna()), "pvalue"].values, method=gesd_mc_method)[1]
                df.loc[(df["pvalue.adj"] < gesd_padj_threshold), "significance"] = "sig"
                pval_threshold = gesd_padj_threshold
        ##### Gesd Test mode end #####

        elif mode == "permutation":
        ##### 2. Permutation mode start #####
            from statsmodels.distributions.empirical_distribution import ECDF

            df["pvalue"] = pd.NA
            df["pvalue.adj"] = pd.NA
            df["significance"] = pd.NA

            print("[",datetime.datetime.now(),"]","Simulate genotype background", flush=True)
            gene_sim_df = simulate_gene(kept_gene_df, sim_fold=fold_times, seed=seed)
            gene_sim_df = gene_sim_df[species_list]
            gene_sim_list = gene_sim_df.index.values
            gene_sim_matrix = gene_sim_df.values

            print("[",datetime.datetime.now(),"]","Perform permutation test", flush=True)
            sim_df = Evolink_calculation(gene_sim_matrix, trait_matrix, species_list, gene_sim_list, tree_file)
            bg = sim_df["Evolink_index"].values
            fg = df.loc[kept_genes, "Evolink_index"].values
            ecdf = ECDF(1-np.abs(bg)) # Empirical cumulative distribution function
            pvals = ecdf(1-np.abs(fg)) # N_bg = len(bg); (ecdf(1-np.abs(fg))*N_bg+1)/(N_bg+1) if to rescale pvalue and make it never be zero
            df.loc[kept_genes, "pvalue"] = pvals
            df.loc[kept_genes, "pvalue.adj"] = smm.multipletests(df.loc[kept_genes, "pvalue"].values, method=perm_mc_method)[1]
            df.loc[(df["pvalue.adj"] < perm_padj_threshold), "significance"] = "sig"
            pval_threshold = perm_padj_threshold
        ##### Permutation mode end #####

        else:
        ##### 3. Z-score mode start #####
            import scipy.stats
            # get Z-score threshold
            # ppf() takes a percentage and returns a standard deviation multiplier for what value that percentage occurs at
            Z_alpha_2 = scipy.stats.norm.ppf(1-alpha/2)
            sd = df.loc[kept_genes, "Evolink_index"].std()
            df["significance"] = pd.NA
            df.loc[kept_genes, "z_score"] = df.loc[kept_genes, "Evolink_index"]/sd
            df.loc[kept_genes, "significance"] = np.where((df.loc[kept_genes, "z_score"].abs() >  Z_alpha_2), "sig", pd.NA)
        ##### Z-score mode end #####

    elif mode == "isolation_forest":
    ##### 4. Isolation Forest mode start #####
        from sklearn.ensemble import IsolationForest
        X = df[["Evolink_index"]]
        
        clf = IsolationForest(n_estimators=n_estimators, max_samples=max_samples)
        _ = clf.fit_predict(X)
        scores = clf.score_samples(X)*(-1)

        ## determine the score_threshold by finding plateau on the curve
        # x = np.arange(0.5, 1, threshold_interval)
        # filtered_scores = scores[scores >= 0.5]
        # outlier_ct = [len([score for score in filtered_scores if score >= i]) for i in x]
        # y = (np.array(outlier_ct) - np.min(outlier_ct)) / (np.max(outlier_ct) - np.min(outlier_ct))
        # slope_Y = np.diff(y)/np.diff(x) # get slopes
        # idx = np.where(slope_Y >= slope_cutoff)[0][0] #-0.01
        # score_thershold = x[idx]
        # index_list = X.index[np.where(filtered_scores >= score_threshold)]

        index_list = X.index[np.where(scores >= outlier_score_threshold)]
        df["significance"] = pd.NA
        df.loc[index_list, "significance"] = "sig"
    ##### Isolation Forest mode end #####

    elif mode == "cutoff":
    ##### 5. Cutoff mode start #####
        thresh = e_threshold
        df["significance"] = pd.NA
        df.loc[df["Evolink_index"].abs() >= thresh, "significance"] = "sig"
    ##### Cutoff mode end #####

    df.reset_index(inplace=True)
    df = df.rename(columns = {'index':'orthoID'})
    res_df = df
    sig_gene_ct = res_df[res_df["significance"]=="sig"].shape[0]
    print("[",datetime.datetime.now(),"]","Find", sig_gene_ct, "significant genes", flush=True)
    
    ### Output ###
    print("[",datetime.datetime.now(),"]","Output result", flush=True)
    outfile = os.path.join(output, "result.tsv")
    res_df.to_csv(outfile, sep="\t", index=False, na_rep='NA')

    ### Change plot according to mode ###
    ### Plot ####
    if plot:
        print("[",datetime.datetime.now(),"]","Generate figures", flush=True)
        iTOL_input(tree_file, trait_df, gene_df, res_df, output, pos_top, neg_top)
        plot_fig(mode, gene_table, trait_table, tree_file, output, pval_threshold, pos_top, neg_top, display_mode)

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Evolink is designed to find gene families associated with trait by explicitly using phylogeny information.')
    
    # Essential Input
    parser.add_argument('-g', '--genotype', help='Tab-delimited gene presence/absence or copy number table. Columns are gene families, while rows are tip names/species/genomes in the phylogenetic tree. If copy number table is provided, please use -c option so that it will be internally converted to binary table. Presence=1, Absence=0.', required=True, dest='gene_table', metavar='GENE_TABLE')
    parser.add_argument('-t', '--phenotype', help='Two-column (so far only one trait is allowed each time) tab-delimited trait presence/absence table. The first column is tip names and the second column is the presence/absence of this trait on the tips/species/genomes. Presence=1, Absence=0.',required=True, dest='trait_table', metavar='TRAIT_TABLE')
    parser.add_argument('-n', '--phylogeny', help='A phylogentic tree in newick format. The tip names should be the same in the gene table and trait table.',required=True, dest='tree_file', metavar='TREE')

    # Optional Input
    parser.add_argument('-m', '--mode', help='Evolink has 4 modes to detect siginificant genotypes assoicated with phenotype: gesd_test, isolation_forest, permutation, z_score and cutoff. Running time: permutation > isolation_forest > gesd_test > z_score > cutoff', required=False, default='gesd_test', dest='mode', choices=["gesd_test", "isolation_forest", "permutation", "z_score", "cutoff"] , metavar='MODE')
    parser.add_argument('-c', '--copy-number', help='The given gene table stores numbers (e.g. gene copy numbers) instead of presence/absence binary values. [Default: True]', action='store_true', required=False, default=False, dest='CN')
    parser.add_argument('-p', '--p-threshold', help='Absolute Prevalence index threshold to filter genes and get Evolink index distribution [Range: 0-1; Default: 0.9]', required=False, default=0.9, type=float, dest='p_threshold', metavar='THRESHOLD')
    parser.add_argument('-r', '--seed', help='Set seed for simulation for reproducibility of the results [Default: 1]', required=False, default=1, type=int, dest='seed', metavar='SEED')
    
    # Gesd Test Mode Options
    parser.add_argument('--gesd-mc-method', help='Multitest correction method [Choices: none, bonferroni, fdr, holm, hommel; Default: none]', choices = ["none", "fdr_bh", "bonferroni", "holm", "hommel"], required=False, type=str, default="none", dest='gesd_mc_method', metavar='METHOD')
    parser.add_argument('--gesd-pval-threshold', help='Original p-value threshold [Default:0.1]', required=False, type=float, default=0.1, dest='gesd_pval_threshold', metavar='THRESHOLD')
    parser.add_argument('--gesd-padj-threshold', help='Adjusted p-value threshold [Default:0.2]', required=False, type=float, default=0.2, dest='gesd_padj_threshold', metavar='THRESHOLD')

    # Cutoff Mode Option
    parser.add_argument('-e', '--e-threshold', help='Absolute Evolink index threshold to select significant genes. [Range: 0-1; Default: 0.375]', required=False, default=0.375, type=float, dest='e_threshold', metavar='THRESHOLD')

    # Isolation Forest Mode Options
    parser.add_argument('--outlier-score-threshold', help='A threshold to determine outliers by IsolationForest [Default: 0.8; Range: 0.5-1]', required=False, default=0.8, type=float, dest='outlier_score_threshold', metavar='THRESHOLD')
    parser.add_argument('--n-estimators', help='Number of tree estimators used in IsolationForest [Default: 200]', required=False, type=int, default=200, dest='n_estimators', metavar='NUMBER')
    parser.add_argument('--max-samples', help='Percentage of training samples for each tree in IsolationForest [Default: 0.1]', required=False, type=float, default=0.1, dest='max_samples', metavar='PERCENTAGE')
    # parser.add_argument('--slope-cutoff', help='A slope cutoff to find the plateau of the IsolationForest score curve [Default: -0.001]', required=False, default=-0.01, type=float, dest='slope_cutoff', metavar='THRESHOLD')
    # parser.add_argument('--threshold-interval', help='The interval/resolution to find the best threshold in IsolationForest scores [Default: 0.005]', required=False, default=0.005, type=float, dest='threshold_interval', metavar='INTERVAL')
    
    # Permutation Mode Options
    parser.add_argument('--fold-times', help='Simulate N*the bumber of genotype input provided by users after filtering (namely simulate N*nrow(gene matrix after filtering) times) [Default: 10]', required=False, default=10, type=int, dest='fold_times', metavar='FOLD_TIMES')
    parser.add_argument('--perm-mc-method', help='Multitest correction method [Choices: bonferroni, fdr_bh, holm, hommel; Default: fdr_bh]', choices = ["fdr_bh", "bonferroni", "holm", "hommel"], required=False, type=str, default="fdr_bh", dest='perm_mc_method', metavar='METHOD')
    parser.add_argument('--perm-padj-threshold', help='Adjusted p-value threshold [Default:0.001]', required=False, type=float, default=0.001, dest='perm_padj_threshold', metavar='THRESHOLD')
    # parser.add_argument('-@', '--threads', help='Threads for permutation test [Default: 4]', required=False, default=4, type=int, dest='threads', metavar='THREADS')

    # Z-score Mode Option => modified z-score?
    parser.add_argument('-a', '--alpha', help='Alpha threshold [Default:0.001; Range: 0.1-0.0001]', required=False, type=float, default=0.001, dest='alpha', metavar='ALPHA')
    
    # Plot Options
    parser.add_argument('-v', '--visualization', help='Whether to generate plots', action='store_true', required=False, default=False, dest='plot')
    parser.add_argument('-N', '--top-genes', help='Top positively and negatively associated genes mapped to tree. [Default: 5,5 for top 5 pos genes and top 5 neg genes.]', required=False, default='5,5', type=str, dest='top_genes')
    parser.add_argument('-d', '--display-mode', help='Tree display mode. [1: circular, 2: rectangular; Default: 1]', type=int, choices=[1, 2], required=False, default=1, dest='display_mode')
    parser.add_argument('-f', '--force', help='Force to overwrite output folder. [Default: False]', action='store_true', required=False, default=False, dest='force')
    
    # Output
    parser.add_argument('-o', '--output', help='output directory', required=True, dest='output', metavar='OUTPUT')
    
    args = parser.parse_args()
    pipeline(args)
