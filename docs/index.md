# Welcome to Evolink

## Overview
In breif, Evolink is a phylogeny-based tool to detect genes (both positively and negatively associated ones) contributed to a phenotype present in multi-species (e.g. resistance, virulence, host and colony).

Identification of genotyep-phenotype associations is a fundamental task not only in microbiology but also in the whole field of biology. Yet as microbial data is rapidly increasing, the scales of gene family pool (~10^6) and phylogenetic tree (with > 10^5 leaves) make current methods less efficient to link genes to traits. 

Phylogenetic information is accepted as a good resource to control for population structure in microbial genotyep-phenotype association analyses, avoiding spurious findings. That's why Evolink was developed based on the use of phylogeny.

Tested on a self-made flagella dataset with a large tree (with 1,948 leaves) and a gene family presence/absence matrix (containing 149,316 gene families), Evolink could give results in less than 5 minutes, demonstrating its capability of mining genes correlated to a phenotype on large-scale datasets.

## Requirement
- Install [Anaconda](https://www.anaconda.com/products/distribution) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) -- A distribution of the Python and R programming languages for scientific computing, greatly simplifying package management and deployment.
- [Anaconda or Miniconda?](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html#anaconda-or-miniconda)

## Installation
To install Evolink is easy.
- Step 1. Git clone project
`git clone https://github.com/nlm-irp-jianglab/Evolink.git`  
`cd Evolink`

- Step 2. Build conda environment
`conda env create -f environment.yml`

- Step 3. Activate Evolink environment
`conda activate Evolink`

- Step 4. Setup R packages
`Rscript setup.R`

## Input
Evolink takes 3 essential input files:
1. Species tree (newick format). It is recommneded that the tree is rooted. Internal node names could be omitted.
> e.g. (species_1:1,(species_2:1,(species_3:1,species_4:1)Internal_1:0.5)Internal_2:0.5)Root:0.1;

2. Trait/Phenotype binary file (tab separated file). The header is a must and should be "Tip" and "Status". Tip column contains the tip names the same as the tree, while Status column contains the presence (1) and absence (0) of the phenotype for each leaf. So far only 1 or 0 is accepted and all leaves should be labeled with a 0/1 status. For example:  

| Tip       | Status |
|-----------|--------|
| species_1 | 0      |
| species_2 | 1      |
| species_3 | 1      |
| species_4 | 0      |

1. Gene presence/absence matrix file (tab separated file). Each row is the binary (0/1) status of each gene cross all species. Each gene should appear in a species for at least one time. The first colname could be any word, but "orthoID" (orthogroup ID) is a nice choice to be shown here. For example:  

| orthoID | species_1 | species_2 | species_3 | species_4 |
|---------|-----------|-----------|-----------|-----------|
| gene_1  | 0         | 1         | 1         | 0         |
| gene_2  | 1         | 1         | 0         | 0         |
| gene_3  | 1         | 0         | 1         | 0         |
| gene_4  | 0         | 1         | 0         | 1         |

## Usage
```
usage: Evolink.py [-h] -g GENE_TABLE -t TRAIT_TABLE -n TREE [-c] 
[-p THRESHOLD] [-e THRESHOLD] [-a ALPHA] [-v] [-N TOP_GENES] [-m {1,2}] [-f] -o OUTPUT

Evolink is designed to find gene families associated with given trait 
with the help ofphylogeny information.

optional arguments:
  -h, --help            show this help message and exit
  -g GENE_TABLE, --genotype GENE_TABLE
                        Tab-delimited gene presence/absence or copy number table. 
                        Columns are gene families, while rows are tip 
                        names/species/genomes in the phylogenetic tree. 
                        If copy number table is provided, please use -c option 
                        so that it will be internally converted to binary table. 
                        Presence=1, Absence=0.
  -t TRAIT_TABLE, --phenotype TRAIT_TABLE
                        Two-column (so far only one trait is allowed each time) 
                        tab-delimited trait presence/absence table. 
                        The first column is tip names and the second column is 
                        the presence/absence of this trait on the 
                        tips/species/genomes. Presence=1, Absence=0.
  -n TREE, --phylogeny TREE
                        A phylogentic tree in newick format. 
                        The tip names should be the same in the gene table and 
                        trait table.
  -c, --copy_number     The given gene table stores numbers (e.g. gene copy numbers) 
                        instead of presence/absence binary values. 
                        [Default: True]
  -p THRESHOLD, --p_threshold THRESHOLD
                        Absolute Prevalence index threshold to filter in genes 
                        for permutation tests [Range: 0-1; Default: 0.9]
  -e THRESHOLD, --e_threshold THRESHOLD
                        Absolute Evolink index threshold to filter in genes for 
                        permutation tests [Range: 0-1; Default: 0.1]
  -a ALPHA, --alpha ALPHA
                        Tail of area cutoff 
                        [Range: 0.8, 0.85, 0.90, 0.95, 0.99, 0.995, 0.997, 0.999; 
                        Default: 0.997]
  -v, --visualization   Whether to generate plots
  -N TOP_GENES, --top_genes TOP_GENES
                        Top positively and negatively associated genes mapped to tree. 
                        [Default: 5,5 for top 5 pos genes and top 5 neg genes.]
  -m {1,2}, --display-mode {1,2}
                        Tree display mode. [1: circular, 2: rectangular; Default: 1]
  -f, --force           Force to overwrite output folder. [Default: False]
  -o OUTPUT, --output OUTPUT
                        output directory
```

## Examples
1. With binary gene presence/absence table (most common usage) and no plots by default:
```
python Evolink.py -g test/gene.tsv -t test/trait.tsv -n test/tree.nwk -o output_dir
```

2. With gene copy number table (add "-c" option):
```
python Evolink.py -g test/gene_CN.tsv -c -t test/trait.tsv -n test/tree.nwk -o output_dir_CN
```

3. Enable plot function (add "-v" option. To save time, Evolink will not generate figures by default):
```
python Evolink.py -g test/gene.tsv -c -t test/trait.tsv -n test/tree.nwk -o output_dir -v
```

4. More complex usage. "-N" is to map top nine positively and top eight negatively associated genes in the plot; "-m" is to use circular layout for the tree and "-f" is to force ovrewrite the output directory if it already exists:
```
python Evolink.py -g test/gene.tsv -c -t test/trait.tsv -n test/tree.nwk -o output_dir -v -N 9,8 -m 1 -f
```

## Output
- A basic output file from Evolink is named "result.tsv" in the output directory provided by the user. It includes "Evolink_index", "Prevelance_index", "significance" and "z_score".  "Evolink_index" and "significance" are the most useful values. For example:  

| orthoID | Prevalence_index     | Evolink_index        | significance | z_score             |
|---------|----------------------|----------------------|--------------|---------------------|
| COG1797 | 0.4759221580503997   | 0.5240778359135583   | sig          | 3.0518889144100863  |
| COG0017 | -0.21711438404978406 | -0.2068655103280859  | NA           | -1.2046503677140863 |
| COG0574 | 0.9447156364996409   | 0.05528435746431705  | NA           | 0.32194018926887563 |
| COG2359 | -0.24181029630947368 | -0.5551479548794538  | sig          | -3.232821106431004  |
| COG2508 | 0.25666479445273405  | 0.530895338741524    | sig          | 3.091589622737808   |
| COG1804 | 0.21892456915043446  | 0.28915028879254645  | NA           | 1.6838234714241953  |
| COG1344 | 0.06105077581463875  | -0.16866160645470016 | NA           | -0.9821756459675905 |
| COG1977 | -0.18435602669012327 | -0.6126022244988042  | sig          | -3.567397454677384  |

- When enabling the plot function (with -v or --visualization option), Evolink provides in the output directory four types of figures:
### iTOL website input
a tree file (input.tree) and annotation file (binary.txt) as well as a zipped file called **Evolink_itol_input.zip** for users to visualize their trees on the [Tree of Life (iTOL)](https://itol.embl.de/).

1. You can simply upload the **input.tree** to iTOL website and drag **binary.txt** to the tree page for visualization.
2. Or after installing [iTOL API](https://github.com/iBiology/iTOL), you can use the following command line to upload and annotate your tree if you have a subscription API key, considering iTOL is not free for batch upload:
`itol Evolink_itol_input.zip -i <your iTOL upload API key> -p <project_name>`

In addition, we also provided a script "Evolink_plot.R" to individually generate local figures after you get the result from Evolink (i.e. result.tsv):
`Rscript --vanilla Evolink_plot.R -g test/gene.tsv -t test/trait.tsv -n test/tree.nwk -r test/output_dir/result.tsv -o test/plot_dir`

### ggtree plot for positively and negatively associated genes
**Tree_mapping_pos.pdf**
![ggtree_pos](images/ggtree_pos.png)

**Tree_mapping_neg.pdf**
![ggtree_neg](images/ggtree_neg.png)

### Evolink plot
**Evolink.pdf**
![Evolink_plot](images/Evolink_plot.png)

### Manhattan plot
**Manhattan.pdf**
![Manhattan_plot](images/Manhattan_plot.png)
