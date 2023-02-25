library(tidyverse)
library(psych)
library(foreach)
library(PRROC)
library(tidymodels)
library(doParallel)
registerDoParallel(32)
set.seed(1234567)
# Usage: time Rscript --vanilla run_tetcorr.R gene.tsv trait.tsv output

######################################################
# Wrap Tetrachoric correlation running in a R script #
######################################################

args = commandArgs(trailingOnly=TRUE)

### read data ###
gene_path = args[1] # row: gene, col: tip
phen_path = args[2]
outfile = args[3] # output dir

gene = read_tsv(gene_path, show_col_types = FALSE)
trait = read_tsv(phen_path, show_col_types = FALSE)
col1_name = colnames(gene)[1]
gene <- gene %>% select(all_of(col1_name), trait$Tip) # reorder gene df

# running test
start=1
end=nrow(gene)
res <- foreach(i=start:end, .combine=rbind) %dopar% {
    orthoID = gene$orthoID[i]
    # print(paste0(i, ":", orthoID))
    row=as.numeric(gene[i,-1])
    if(sum(row==1)==length(row) | sum(row==0)==length(row)){corr = NA}
    else{
        dt = data.frame(gene=row, trait=trait$Status)
        corr = tetrachoric(dt, correct=0)$rho[1,2]
    }
    c(orthoID, corr)
}

res <- res %>% as.data.frame()
colnames(res) <- c("orthoID", "Tetrachoric_correlation")
res$Tetrachoric_correlation = as.numeric(res$Tetrachoric_correlation)
# res$Pvalue = as.numeric(res$Pvalue)
# res$Padj = p.adjust(res$Pvalue, method="fdr")
write_tsv(res, paste0(args[3], "/output.txt"))

############### get optimal cutoff using ground truth ##############
gold_genes=c(paste0("neg", 1:10), paste0("pos", 1:10))
PR_dat = function(dt, gold_genes, colname){
    dt$score = abs(dt[[colname]])
    dt = dt %>% mutate(true_class=ifelse(orthoID %in% gold_genes, "positive", "negative")) %>%
                     column_to_rownames(var="orthoID") %>%
                     select(true_class, score)
    dt$true_class = factor(dt$true_class, levels=c("positive", "negative"))
    pr_dt <- dt %>% pr_curve(true_class, score)
    return(pr_dt)
}
rownames(res)=NULL
pr_dt = PR_dat(res, gold_genes, "Tetrachoric_correlation")
pr_dt=pr_dt %>% mutate(F1=(2*precision*recall)/(precision + recall)) %>% arrange(desc(F1))
cutoff = pr_dt[1,][[".threshold"]]
print(paste0("correlation cutoff:", cutoff))

res = res %>% mutate(sig=ifelse((abs(Tetrachoric_correlation) >= cutoff) & !is.na(Tetrachoric_correlation), "sig", "background"))
    tetcorr_sig_genes = res %>% filter(sig=="sig") %>% select(orthoID)
    write_tsv(tetcorr_sig_genes, paste0(args[3], "/Tetcorr_sig_genes.list"), col_names=FALSE)
######################################################