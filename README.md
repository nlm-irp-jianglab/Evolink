# Evolink

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

**Evolink** is a phylogenetic approach for rapid identification of genotype-phenotype associations in large-scale microbial multi-species data.

![Evolink](img/Logo.jpg)

# File structure and script description
```
.
├── README.md
├── Evolink.py
├── Evolink_plot.R
├── scripts
│   ├── gesd_test.R
│   └── simulate_gene.R
├── environment.yml
├── setup.R
├── docs/
├── img/
├── flagella_data/
├── gram_staining_data/
├── test_data/
└── Evolink_paper/
```

- `README.md`: The page you are seeing now.
- `Evolink.py`: Evolink main script.
- `Evolink_plot.R`: An indenpendent R script to help visualize Evolink results.
- `scripts/gesd_test.R`: An R script to perform the gesd test to detect phenotype-associated genotypes (Evolink mode=gesd_test).
- `scripts/simulate_gene.R`: An R script to perform permutation on genotypes to detect phenotype-associated genotypes (Evolink mode=permutation).
- `environment.yml`: Conda environment yaml file with all dependencies for creating an Evolink environment. See details in [Evolink document](https://nlm-irp-jianglab.github.io/Evolink).
- `setup.R`: An R script to help install R packages used in Evolink. See details in [Evolink document](https://nlm-irp-jianglab.github.io/Evolink).
- `docs/`: A folder containing materials for [Evolink document](https://nlm-irp-jianglab.github.io/Evolink).
- `img/`: A folder containing the Evolink logo.
- `flagella_data/`: A folder containing a walkthrough named `Records.md` and related data for a use case with flagella data as an example. See details in [Flagella use case](https://github.com/nlm-irp-jianglab/Evolink/blob/main/flagella_data/Records.md).
- `gram_staining_data/`: A folder containing a walkthrough named `Records.md` and related data for a use case with gram-staining data as an example.
- `test_data/`: A folder containing a bash script `test.sh` and related data for testing.
- `Evolink_paper/`: A folder containing additional data used in the Evolink paper along with scripts to compare Evolink to other methods on simulated and empirical datasets. See details in `Evolink_paper/README.md`.

# Document

Documentation detailing the installation, Evolink input and Evolink output is provided here: [Evolink document](https://nlm-irp-jianglab.github.io/Evolink).

# Flagellar dataset example

A use case showing how to apply Evolink to a real-world flagella-related dataset is provided here: [Flagella use case](https://github.com/nlm-irp-jianglab/Evolink/blob/main/flagella_data/Records.md). The genome, gff annotation and emapper annotation files are provided at 

# Gram-staining dataset example

A use case showing users how to apply Evolink to a real-world gram-staining dataset is provided here: [Gram-staining use case](https://github.com/nlm-irp-jianglab/Evolink/blob/main/gram_staining_data/Records.md).

# Other resources

- We also provide a [docker](https://hub.docker.com/r/nlmirpjianglab/evolink) and a [web portal](https://jianglabnlm.com/evolink) to make it easier to run and access Evolink on different platforms.  

- The genome, gff annotation and emapper annotation files for Flagellar and Gram-staining datasets are provided at [WOL_flagellum_data.tar.gz](https://ftp.ncbi.nlm.nih.gov/pub/mgx/WOL/WOL_flagellum_data.tar.gz) and [WOL_gram_staining_data.tar.gz](https://ftp.ncbi.nlm.nih.gov/pub/mgx/WOL/WOL_gram_staining_data.tar.gz), respectively.

# Acknowledgements

- [Unifrac](https://github.com/biocore/unifrac) python module
- [Scikit-learn](https://scikit-learn.org/stable/) python library
- [PMCMRplus](https://cran.r-project.org/web/packages/PMCMRplus/index.html) R package
- [treeWAS](https://github.com/caitiecollins/treeWAS) R package

# Citation
Yiyan Yang , Xiaofang Jiang, Evolink: a phylogenetic approach for rapid identification of genotype–phenotype associations in large-scale microbial multispecies data, Bioinformatics, Volume 39, Issue 5, May 2023, btad215, https://doi.org/10.1093/bioinformatics/btad215
