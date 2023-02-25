# Data and scripts for Evolink paper
This a pipeline for generating simulated data and runing Evolink and alternative methods on both similated data and empirical data.

## Install environment
1. git clone project  
`git clone https://github.com/nlm-irp-jianglab/Evolink.git`  
`cd Evolink/Evolink_paper`  

2. install conda env  
`conda env create -f env_yml/bio-env.yml`  
`conda env create -f env_yml/pyseer_env.yml`  

3. install R packages  
`conda activate bio-env`  
Type `R` and enter R terminal:  
`install.packages(c("tidyverse", "phytools", "geiger", "snow", "VGAM", "foreach", "doParallel"))`  

Note: We highly recommend users who want to follow this protocol run scripts on an HPC (high performance cluster) node with the [Slurm Workload Manager](https://slurm.schedmd.com/overview.html) on Linux system. Because the memory needed to generate simulated data will exceed 600 GB with a maximal memery of 1500 GB, considering that most of the simulated datasets have >= 1000 species and >= 20K genes.

## Run pipeline
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
After running, each method will have a `<method_name>_out` folder in each simulated dataset folder containing all the outputs of this method. One of the most important outputs is the file `<method>_sig_genes.list` listing the significant gene families found by this method.

### Step.3 Run tested methods on flagella data
```
/bin/bash Step3.run_methods_on_flagella.sh
```
After running, each method will have a `<method_name>_out` folder in the flagella (sub)dataset folder containing all the outputs of this method and the `<method>_sig_genes.list` file within the folder lists the significant gene families found by this method.

### Step.4 Run tested methods on gram-staining data
```
/bin/bash Step4.run_methods_on_gram-staining.sh
```
After running, each method will have a `<method_name>_out` folder in the gram stainig (sub)dataset folder containing all the outputs of this method and the `<method>_sig_genes.list` file within the folder lists the significant gene families found by this method.
