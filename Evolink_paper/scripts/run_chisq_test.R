library(tidyverse)

# Usage: time Rscript --vanilla run_chisq_test.R gene.tsv trait.tsv

my_chisq <- function(y){
  # include here as.numeric to be sure that your values are numeric:
  table <- matrix(as.numeric(c(y[2], y[3], y[4], y[5])), ncol = 2, byrow = TRUE)
  p <- chisq.test(table, correct=F)
  return(list(x=p$p.value, y=p$statistic))
}

args = commandArgs(trailingOnly=TRUE)

### read data ###
gene_path = args[1]
phen_path = args[2]

gene = read_tsv(gene_path, show_col_types = FALSE)
trait = read_tsv(phen_path, show_col_types = FALSE)
col1_name = colnames(gene)[1]
gene <- gene %>% select(all_of(col1_name), trait$Tip) # reorder gene df

gene_mat <- gene[,-1]
trait_mat <- trait[,-1]
A = sweep(gene_mat, MARGIN=2, t(trait_mat), `*`) # G1P1
B = sweep(1-gene_mat, MARGIN=2, t(trait_mat), `*`) # G0P1
C = sweep(gene_mat, MARGIN=2, t(1-trait_mat), `*`) # G1P0
D = sweep(1-gene_mat, MARGIN=2, t(1-trait_mat), `*`) # G0P0

a = A %>% mutate(rowsum = rowSums(.))
b = B %>% mutate(rowsum = rowSums(.))
c = C %>% mutate(rowsum = rowSums(.))
d = D %>% mutate(rowsum = rowSums(.))

df = cbind.data.frame(orthoID=gene[[col1_name]], G1P1=a$rowsum, G0P1=b$rowsum, G1P0=c$rowsum, G0P0=d$rowsum)

# running test
dt = df
options(warn = -1)
res <- apply(dt, 1, my_chisq)
defaultW <- getOption("warn")
dt$chisq.pval=unlist(res)[grepl('x', names(unlist(res)))]
dt$chisq.Xsquared=unlist(res)[grepl('X-squared', names(unlist(res)))]
options(warn = defaultW)

dt1 = dt
dt1$chisq.padj = p.adjust(dt1$chisq.pval, method="fdr")
dt1 = dt1 %>% mutate(sig=ifelse((chisq.padj < 0.05), "sig", "background"))
chisq_sig_genes = dt1 %>% filter(sig=="sig") %>% select(orthoID)
write_tsv(chisq_sig_genes, args[3], col_names=FALSE)