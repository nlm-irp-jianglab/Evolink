args = commandArgs(trailingOnly=TRUE)
outdir = args[1]

load(paste0(outdir, "/hogwash_phyc_", outdir, ".rda"))
load(paste0(outdir, "/hogwash_synchronous_", outdir, ".rda"))
x=hogwash_phyc$sig_pvals
y=hogwash_synchronous$sig_pvals
write.table(row.names(x), paste0(outdir, "/hogwash_phyc_sig_genes.list"), quote=F, col.names=F, row.names=F)
write.table(row.names(y), paste0(outdir, "/hogwash_synchronous_sig_genes.list"), quote=F, col.names=F, row.names=F)