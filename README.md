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
└── test_data/

```
- `README.md`: The page you are watching now.
- `Evolink.py`: Evolink main script.
- `Evolink_plot.R`: An indenpendent R script to help visualize Evolink result.
- `scripts/gesd_test.R`: An R script to perform the gesd test to detect phenotype-associated genotypes (Evolink mode=gesd_test).
- `scripts/simulate_gene.R`: An R script to perform permutation on genotype to detect phenotype-associated genotypes (Evolink mode=permutation).
- `environment.yml`: Conda environment yaml file with all dependencies for creating Evolink environment. See details in [Evolink document](https://nlm-irp-jianglab.github.io/Evolink).
- `setup.R`: An R script to help install R packages used in Evolink. See details in [Evolink document](https://nlm-irp-jianglab.github.io/Evolink).
- `docs/`: A folder containing materials for [Evolink document](https://nlm-irp-jianglab.github.io/Evolink).
- `img/`: A folder containing the Evolink logo.
- `flagella_data/`: A folder containing a walkthrough named `Records.md` and related data for a use case with flagella data as an example. See details in [Flagella use case](https://github.com/nlm-irp-jianglab/Evolink/blob/main/flagella_data/Records.md).
- `gram_staining_data/`: A folder containing a walkthrough named `Records.md` and related data for a use case with gram-staining data as an example.
- `test_data/`: A folder containing a bash script `test.sh` and related data for testing.


# Document

A detailed guidance is provided here: [Evolink document](https://nlm-irp-jianglab.github.io/Evolink).

# A use case with flagella data as an example

A use case showing users how to apply Evolink to a real-world flagella-related dataset is provided here: [Flagella use case](https://github.com/nlm-irp-jianglab/Evolink/blob/main/flagella_data/Records.md).

# Other resources

We also provide a [docker](https://hub.docker.com/r/nlmirpjianglab/evolink) and a [web portal](https://jianglabnlm.com/evolink) to make you have a better experience with Evolink.

# Acknowledgements

- [Unifrac](https://github.com/biocore/unifrac) python module
- [Scikit-learn](https://scikit-learn.org/stable/) python library
- [PMCMRplus](https://cran.r-project.org/web/packages/PMCMRplus/index.html) R package
- [treeWAS](https://github.com/caitiecollins/treeWAS) R package

# Citation
To be added
