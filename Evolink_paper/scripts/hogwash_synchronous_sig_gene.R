args = commandArgs(trailingOnly=TRUE)
outdir = args[1]

load(paste0(outdir, "/hogwash_synchronous_", outdir, ".rda"))
y=hogwash_synchronous$sig_pvals
write.table(row.names(y), paste0(outdir, "/hogwash_synchronous_sig_genes.list"), quote=F, col.names=F, row.names=F)
z=hogwash_synchronous$raw_pvals
write.table(data.frame(orthoID=rownames(z), raw_pval=z), paste0(outdir, "/hogwash_synchronous_rawpval_genes.list"), quote=F, row.names=F)
