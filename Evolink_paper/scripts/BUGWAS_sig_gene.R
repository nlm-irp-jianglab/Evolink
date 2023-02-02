args = commandArgs(trailingOnly=TRUE)
outdir = args[1]

dt = read.table(paste0(outdir, "/bugwas_output_biallelic_lmmout_allSNPs.txt"), header=T)
dt$qvalue = p.adjust(dt$p_lrt, method="fdr")
dt1 = dt[(dt$qvalue<0.05 & !is.na(dt$qvalue) & dt$qvalue>=0), "ps"]
write.table(dt1, paste0(outdir, "/BUGWAS_sig_genes.list"), quote=F, col.names=F, row.names=F)
