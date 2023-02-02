suppressMessages({
    library(tidyverse)
    library(treeWAS)
})

# Version v1.0: filter genes that are too prevalent and rare, so as to reduce running time
# Usage: time Rscript --vanilla run_treeWAS.R tree.nwk gene.tsv trait.tsv treeWAS_out

args = commandArgs(trailingOnly=TRUE)

### read data ###
tree_path = args[1]
gene_path = args[2]
phen_path = args[3]
outdir = args[4]

if(!file.exists(outdir)){dir.create(outdir)}

phylo <- read.tree(tree_path) # need internal nodes if provide ASR inputs
gene <- read.table(gene_path, header=T, row.names=1, comment.char="", sep="\t") # rows: genes; cols: tips, first colname=tip
phen <- read.table(phen_path, header=T) # header: Tip, Status

### reorder dataframe by tree ###
tip_names <- phylo$tip.label
# node_names <- c(phylo$tip.label, phylo$node.label)

genePA = as.data.frame(t(gene)) # transpose gene dataframe, cols: genes; rows: tips
genePA1 <- genePA[order(match(rownames(genePA), tip_names)),] # reorder genePA rows
phen1 <- phen[order(match(phen$Tip, tip_names)),]
phen2 <- factor(phen1$Status, levels=c("0","1"))
names(phen2) <- phen1$Tip

### run treeWAS ###
out <- treeWAS(snps = genePA1,
                phen = phen2,
                tree = phylo,
                p.value=0.05,
                seed = 1,
                plot.tree = FALSE,
                plot.manhattan = FALSE,
                plot.null.dist = FALSE)

dt = data.frame(orthoID=names(out$terminal$corr.dat), 
              terminal_score=out$terminal$corr.dat,
              terminal_pval=out$terminal$p.vals,
              simultaneous_score=out$simultaneous$corr.dat,
              simultaneous_pval=out$simultaneous$p.vals,
              subsequent_score=out$subsequent$corr.dat,
              subsequent_pval=out$subsequent$p.vals)

df = data.frame(sig_gene=out$treeWAS.combined$treeWAS.combined)
write_tsv(df, paste0(outdir, "/treeWAS_sig_genes.list"), col_names=FALSE)
write.table(dt, file=paste0(outdir, "/treeWAS_output.tsv"), sep="\t", quote=F, row.names=F)
