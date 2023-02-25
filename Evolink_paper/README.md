# Data and scripts for Evolink paper
This a pipeline for generating simulated data and runing Evolink and alternative methods on both similated data and empirical data.

## Prerequisite
Please make sure the following R packages have been installed:
- tidyverse
- phytools
- geiger
- snow
- VGAM
- foreach
- doParallel

We highly recommend users who want to follow this protocol run scripts on an HPC (high performance cluster) node with the Slurm Workload Manager on Linux system. Because the memory needed to generate simulated data will exceed 600 GB with a maximal memery of 1500 GB, considering that most of the simulated datasets have >= 1000 species and >= 20K genes.

## Pipeline
We prepared four bash scripts each representing a step in this pipeline.
```
├── Step1.simulate_data.sh
├── Step2.run_methods_on_simData.sh
├── Step3.run_methods_on_flagella.sh
└── Step4.run_methods_on_gram-staining.sh
```
Running the pipeline is easy by just going through each script:

### Step.1 Generate simulated data
```
/bin/bash Step1.simulate_data.sh
```
We are aware that the computational resource for generating these datasets is difficult to obtain. Therefore, the resulting datasets are zipped and provided in this folder:
```
DataS1_simData_phenotype_prevalence.zip
DataS2_simData_phenotype_phylo_overdispersion.zip
DataS3_simData_species_number.zip
DataS4_simData_gene_number.zip
```

### Step.2 Run tested methods on simulated data
```
/bin/bash Step2.run_methods_on_simData.sh
```

### Step.3 Run tested methods on flagella data
```
/bin/bash Step3.run_methods_on_flagella.sh
```

### Step.4 Run tested methods on gram-staining data
```
/bin/bash Step4.run_methods_on_gram-staining.sh
```
