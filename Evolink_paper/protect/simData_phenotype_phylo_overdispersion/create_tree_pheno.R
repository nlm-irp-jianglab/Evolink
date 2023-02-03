library(phytools)

Ntips = 1000

while(TRUE){
    pheno_prev_list = c()
    tree = pbtree(n=Ntips)
    tree = makeNodeLabel(tree, method = "number", prefix = "int") # assign intnode names

    for(cluster in c(2,4,8,16,32)){
        # fix phenotypes
        getTips = function(tree, node){
            desc = getDescendants(tree, node)
            tips = desc[desc<=length(tree$tip.labels)]
            return(tips)
        }

        hc = as.hclust.phylo(tree)
        tip_clst = cutree(hc, k=cluster)
        pheno = tip_clst %% 2
        pheno_perv = sum(pheno)/length(pheno)
        pheno_prev_list = c(pheno_prev_list, pheno_perv)

        # outfile
        outdir=file.path("simData_phenotype_phylo_overdispersion", paste0("cluster_", cluster))
        if(!file.exists(outdir)){dir.create(outdir)}

        saveRDS(tree, paste0(outdir, "/tree.rds"))
        saveRDS(pheno, paste0(outdir, "/pheno.rds"))
    }
    # make sure phenotype prevalence is not very imbalanced
    if(all(pheno_prev_list < 0.6) && all(pheno_prev_list > 0.4)) break
}
print(pheno_prev_list)