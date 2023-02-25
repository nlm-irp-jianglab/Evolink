#################################################################
# Extract significant gene families from ForwardGenomics result #
#################################################################

args = commandArgs(trailingOnly=TRUE)
infile = args[1]
outdir = args[2]

dt = read.table(infile, header=T, sep="")
dt$qvalue = p.adjust(dt$GLS_Pvalue, method="fdr")
dt1 = gsub("^X", "", dt[(dt$qvalue < 0.05 & !is.na(dt$qvalue) & dt$qvalue>=0), "elementID"])
write.table(dt1, paste0(outdir, "/FG_sig_genes.list"), quote=F, col.names=F, row.names=F)
