library(hogwash)
library(ape)
library(tidyverse)

# Usage: time Rscript --vanilla run_hogwash.R tree.nwk gene.tsv trait.tsv output

args = commandArgs(trailingOnly=TRUE)

### read data ###
tree_path = args[1]
gene_path = args[2]
phen_path = args[3]
outprefix = args[4]

phylo <- read.tree(tree_path) # need internal nodes if provide ASR inputs
phylo <- multi2di(phylo)
phylo$edge.length[phylo$edge.length==0]=min(phylo$edge.length[phylo$edge.length!=0])
phylo$node.label=rep("100", length(phylo$node.label))

gene <- read.table(gene_path, header=T, row.names=1, comment.char="", sep="\t") # rows: genes; cols: tips, first colname=tip
phen <- read_tsv(phen_path) # header: Tip, Status

### reorder dataframe by tree ###
tip_names <- phylo$tip.label
node_names <- c(phylo$tip.label, phylo$node.label)
genePA = as.data.frame(t(gene))
genePA1 <- as.matrix(genePA[order(match(rownames(genePA), tip_names)),])
phen1 = phen[order(match(phen$Tip, tip_names)),]
phen1 = phen1 %>% column_to_rownames(var="Tip") %>% as.matrix()

### run hogwash ###
print("Running hogwash")
hogwash(pheno=phen1, 
        geno=genePA1,  
        tree=phylo,      
        file_name = outprefix,
        dir = ".",
        perm = 1000,
        fdr = 0.05,
        test = "synchronous")
