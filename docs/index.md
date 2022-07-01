## Welcome to Evolink

### Overview
In breif, Evolink is a phylogeny-based tool to detect genes (both positively and negatively associated ones) contributed to a phenotype present in multi-species (e.g. resistance, virulence, host and colony).

Identification of genotyep-phenotype associations is a fundamental task not only in microbiology but also in the whole field of biology. Yet as microbial data is rapidly increasing, the scales of gene family pool (~10^6) and phylogenetic tree (with > 10^5 leaves) make current methods less efficient to link genes to traits. 

Phylogenetic information is accepted as a good resource to control for population structure in microbial genotyep-phenotype association analyses, avoiding spurious findings. That's why Evolink was developed based on the use of phylogeny.

Tested on a self-made flagella dataset with a large tree (with 1,948 leaves) and a gene family presence/absence matrix (containing 149,316 gene families), Evolink could give results in less than 5 minutes, demonstrating its capability of mining genes correlated to a phenotype on large-scale datasets.

### Prerequisite
Prepare [Anaconda](https://www.anaconda.com/products/distribution) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html). [Anaconda or Miniconda?](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html#anaconda-or-miniconda)
Although tests were performed on Linux system, it is likely Evolink could be run on Windows and Mac since conda could be run cross-platform.

### Installation
To install Evolink is easy:
```
# Step 1. Git clone project
git clone https://github.com/nlm-irp-jianglab/Evolink.git
cd Evolink
# Step 2. Build conda environment
conda env create -f environment.yml
# Step 3. Activate Evolink environment
conda activate Evolink
# Step 4. Set up R packages
Rscript setup.R
```

### Input


### Output


### Usage
```
usage: Evolink.py [-h] -g GENE_TABLE -t TRAIT_TABLE -n TREE [-c] [-p THRESHOLD] [-e THRESHOLD] [-a ALPHA] [-v] [-N TOP_GENES] [-m {1,2}] [-f] -o OUTPUT

Evolink is designed to find gene families associated with given trait with the help of phylogeny information.

optional arguments:
  -h, --help            show this help message and exit
  -g GENE_TABLE, --genotype GENE_TABLE
                        Tab-delimited gene presence/absence or copy number table. Columns are gene families, while rows are tip names/species/genomes in the phylogenetic tree. If copy number table is provided, please use -c option so that it will be internally converted to binary table. Presence=1, Absence=0.
  -t TRAIT_TABLE, --phenotype TRAIT_TABLE
                        Two-column (so far only one trait is allowed each time) tab-delimited trait presence/absence table. The first column is tip names
                        and the second column is the presence/absence of this trait on the tips/species/genomes. Presence=1, Absence=0.
  -n TREE, --phylogeny TREE
                        A phylogentic tree in newick format. The tip names should be the same in the gene table and trait table.
  -c, --copy_number     The given gene table stores numbers (e.g. gene copy numbers) instead of presence/absence binary values. [Default: True]
  -p THRESHOLD, --p_threshold THRESHOLD
                        Absolute Prevalence index threshold to filter in genes for permutation tests [Range: 0-1; Default: 0.9]
  -e THRESHOLD, --e_threshold THRESHOLD
                        Absolute Evolink index threshold to filter in genes for permutation tests [Range: 0-1; Default: 0.1]
  -a ALPHA, --alpha ALPHA
                        Tail of area cutoff [Range: 0.8, 0.85, 0.90, 0.95, 0.99, 0.995, 0.997, 0.999; Default: 0.997]
  -v, --visualization   Whether to generate plots
  -N TOP_GENES, --top_genes TOP_GENES
                        Top positively and negatively associated genes mapped to tree. [Default: 5,5 for top 5 pos genes and top 5 neg genes.]
  -m {1,2}, --display-mode {1,2}
                        Tree display mode. [1: circular, 2: rectangular; Default: 1]
  -f, --force           Force to overwrite output folder. [Default: False]
  -o OUTPUT, --output OUTPUT
                        output directory
```

Examples
```
1. With binary gene presence/absence table (most common):
python Evolink.py -g test/gene.tsv -t test/trait.tsv -n test/tree.nwk -o output_dir
```

```
2. With gene copy number table (add "-c" option):
python Evolink.py -g test/gene_CN.tsv -c -t test/trait.tsv -n test/tree.nwk -o output_dir_CN
```

```
3. Enable plot function (add "-v" option. To save time, Evolink will not generate figures by default):
python Evolink.py -g test/gene.tsv -c -t test/trait.tsv -n test/tree.nwk -o output_dir -v
```

```
4. More complex usage. "-N" is to map top 9 positively and top 8 negatively associated genes in the plot; "-m" is to use circular layout for the tree and "-f" is to force ovrewrite the output directory if it already exists:
python Evolink.py -g test/gene.tsv -c -t test/trait.tsv -n test/tree.nwk -o output_dir -v -N 9,8 -m 1 -f
```

### Gallery
When enabling the plot function, Evolink provides in the ourtput directory a tree file (input.tree) and annotation file (binary.txt) as well as a zipped file called "Evolink_itol_input.zip" for users to visualize their trees on the [Tree of Life (iTOL)](https://itol.embl.de/).

1. You can simply upload the "input.tree" to iTOL website and drag "binary.txt" to the tree page for visualization.
2. Or after installing [iTOL API](https://github.com/iBiology/iTOL), you can use the following command line to upload and annotate your tree:
```
itol Evolink_itol_input.zip -i <your iTOL upload API key> -p <project_name>
```

In addition, we also provided a script "Evolink_plot.R" to individually generate local figures after you get the result from Evolink (i.e. result.tsv):
```
Rscript --vanilla ../Evolink_plot.R -g gene.tsv -t trait.tsv -n tree.nwk -r output_dir/result.tsv -o plot_dir
```

### Support or Contact

- Bulleted
- Lists

1. Number
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.
