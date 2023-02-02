library(bugwas)
args = commandArgs(trailingOnly=TRUE)
phylo <- args[1]
gen <- args[2]
pheno <- args[3]
prefix <- args[4]
gem.path <- "softwares/bugwas/gemma/gemma.0.93b"
data <- lin_loc(gen = gen, pheno = pheno, phylo = phylo, prefix = prefix, gem.path = gem.path,
                creatingAllPlots = FALSE, runTriTetrallelic = FALSE)
