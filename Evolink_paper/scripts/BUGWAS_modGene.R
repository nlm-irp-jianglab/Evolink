library(tidyverse)
args = commandArgs(trailingOnly=TRUE)

uniq_len=function(x){return(length(unique(unlist(x))))}

gene_file = args[1]
outfile = args[2]

df = read.table(gene_file, row.names=1, header=T)
res = apply(df, 1, uniq_len)
df1 = df[!(rownames(df) %in% names(res[res==1])), ]
df1 = df1 %>% rownames_to_column(var="ps")
write.table(df1, outfile, quote=F, row.names=F, sep="\t")