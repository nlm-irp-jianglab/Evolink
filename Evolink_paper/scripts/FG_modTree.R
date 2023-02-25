#################################################
# Modify tree file:                             #
# 1) convert it into bifuricated tree           #
# 2) add an extreme number to zero branch length#
# 3) assign internal nodes names                #
# so that ForwardGenomics can use as input      #
#################################################

library(ape)
args = commandArgs(trailingOnly=TRUE)
t=read.tree(args[1])
t1 <-multi2di(t)
t1$edge.length[t1$edge.length==0]=min(t1$edge.length[t1$edge.length!=0])

t1$node.label=paste0("A",1:length(t1$node.label))
write.tree(t1, args[2])
