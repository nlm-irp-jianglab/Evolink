library(phytools)

Ntips = 1000

tree = pbtree(n=Ntips)
tree = makeNodeLabel(tree, method = "number", prefix = "int") # assign intnode names

for(prev in c(0.1, 0.3, 0.5, 0.7, 0.9)){
    # fix phenotypes
    pheno = rep(NA, Ntips)
    pos_index = sample(1:Ntips, round(Ntips*prev), replace=F)
    neg_index = setdiff(1:Ntips, pos_index)
    # print(length(pos_index))
    # print(length(neg_index))
    pheno = replace(pheno, pos_index, 1)
    pheno = replace(pheno, neg_index, 0)
    names(pheno) = tree$tip.label

    # outfile
    outdir=file.path("simData_phenotype_prevalence", paste0("prev_", prev))
    if(!file.exists(outdir)){dir.create(outdir)}

    saveRDS(tree, paste0(outdir, "/tree.rds"))
    saveRDS(pheno, paste0(outdir, "/pheno.rds"))
}
