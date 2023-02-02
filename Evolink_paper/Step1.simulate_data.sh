cd Evolink_paper/

# simulated data with low to high phenotype phylogenetic overdispersion
Rscript --vanilla simData_phenotype_phylo_overdispersion/create_tree_pheno.R
sbatch simData_phenotype_phylo_overdispersion/sbatch_jobs.sh 2
sbatch simData_phenotype_phylo_overdispersion/sbatch_jobs.sh 4
sbatch simData_phenotype_phylo_overdispersion/sbatch_jobs.sh 8
sbatch simData_phenotype_phylo_overdispersion/sbatch_jobs.sh 16
sbatch simData_phenotype_phylo_overdispersion/sbatch_jobs.sh 32

# simulated data with low to high phenotype prevalence
Rscript --vanilla simData_phenotype_prevalence/create_tree_pheno.R
sbatch simData_phenotype_prevalence/sbatch_jobs.sh 0.1
sbatch simData_phenotype_prevalence/sbatch_jobs.sh 0.3
sbatch simData_phenotype_prevalence/sbatch_jobs.sh 0.5
sbatch simData_phenotype_prevalence/sbatch_jobs.sh 0.7
sbatch simData_phenotype_prevalence/sbatch_jobs.sh 0.9

# simulated data with varying genotype numbers
sbatch simData_genotype_number/sbatch_jobs.sh 1000 80 10000 data_s1k_g10k
sbatch simData_genotype_number/sbatch_jobs.sh 1000 80 20000 data_s1k_g20k
sbatch simData_genotype_number/sbatch_jobs.sh 1000 80 40000 data_s1k_g40k
sbatch simData_genotype_number/sbatch_jobs.sh 1000 80 80000 data_s1k_g80k
sbatch simData_genotype_number/sbatch_jobs.sh 1000 80 160000 data_s1k_g160k

# simulated data with varying species numbers
sbatch simData_species_number/sbatch_jobs.sh 100 20 10000 data_s0.1k_g10k
sbatch simData_species_number/sbatch_jobs.sh 200 40 10000 data_s0.2k_g10k
sbatch simData_species_number/sbatch_jobs.sh 400 60 10000 data_s0.4k_g10k
sbatch simData_species_number/sbatch_jobs.sh 800 80 10000 data_s0.8k_g10k
sbatch simData_species_number/sbatch_jobs.sh 1600 100 10000 data_s1.6k_g10k
sbatch simData_species_number/sbatch_jobs.sh 3200 120 10000 data_s3.2k_g10k

### The generated results are zipped and named:
# DataS1_simData_phenotype_prevalence.zip
# DataS2_simData_phenotype_phylo_overdispersion.zip
# DataS3_simData_species_number.zip
# DataS4_simData_genotype_number.zip