#########################################################
# Extract significant gene families from Phylolm result #
#########################################################

args = commandArgs(trailingOnly=TRUE)
infile = args[1]
outdir = args[2]

dt = read.table(infile, header=T, sep="")
dt$qvalue = p.adjust(dt$pred.Pvalue, method="fdr")
dt1 = gsub("^X", "", dt[(dt$qvalue<0.05 & !is.na(dt$qvalue) & abs(dt$pred.Estimate) >= 1 & dt$qvalue>=0), "orthoID"])
write.table(dt1, paste0(outdir, "/Phylolm_sig_genes.list"), quote=F, col.names=F, row.names=F)
