args = commandArgs(trailingOnly=TRUE)
outdir = args[1]

load(paste0(outdir, "/hogwash_phyc_", outdir, ".rda"))
x=hogwash_phyc$sig_pvals
write.table(row.names(x), paste0(outdir, "/hogwash_phyc_sig_genes.list"), quote=F, col.names=F, row.names=F)
y=hogwash_phyc$raw_pvals
write.table(data.frame(orthoID=rownames(y), raw_pval=y), paste0(outdir, "/hogwash_phyc_rawpval_genes.list"), quote=F, row.names=F)
