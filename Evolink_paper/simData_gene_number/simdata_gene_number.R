
suppressMessages({
    library(tidyverse)
    library(phytools)
    library(geiger)
    library(snow)
    library(VGAM)    
    library(foreach)
    library(doParallel)}
)
args = commandArgs(trailingOnly=TRUE)

modify_prob <- function(prob, adjust_prob){
    ################################################
    ### This function is used to adjust prob     ###
    ### based on ancestral phenotype probability ###
    ################################################

    # adjust_prob: c(P(pheno=0), P(pheno=1))
    alpha = prob[2]/prob[4] # get ratio of 01/11
    beta = prob[1]/prob[3] # get ratio of 00/10
    x = adjust_prob[1] # pheno==1 probability
    y = 1 - x
    b_1 = alpha*y/(1+alpha)
    d_1 = y/(1+alpha)
    a_1 = beta*x/(1+beta)
    c_1 = x/(1+beta)
    return(c(a_1, b_1, c_1, d_1))
}

# set parameters #
Ntips = as.numeric(args[1]) # 1000
Nsubs = as.numeric(args[2]) # 80, total number of substitution of characters along the whole tree
associate_factor = 15 # control the correlation degree for geno and pheno
all_ncpus = 72 #72
sim_ncpus = 50 #50 # should be dividable by sim_times
sim_times = 10000 #5000, the times of simmap phenotype on ancestral nodes
N_all_geno = as.numeric(args[3]) # 10000 20000 40000 80000 160000
N_pos_geno = 10
N_neg_geno = 10
registerDoParallel(all_ncpus)

# randomly generate a tree
tree = pbtree(n=Ntips, scale=10)
tree = makeNodeLabel(tree, method = "number", prefix = "int") # assign intnode names

# outfile
outdir = args[4]
# dir.create(outdir)

# design a correlated Q matrix
sub_rate=Nsubs/sum(tree$edge.length)

# pos Q matrix
x1=x2=n1=n2=sub_rate
y1=y2=m1=m2=associate_factor*sub_rate
Q_pos = matrix(c(0,x1,x2,0,
                y1,0,0,y2,
                m1,0,0,m2,
                0,n1,n2,0),
                4,4,byrow=TRUE)
rownames(Q_pos)<-colnames(Q_pos)<-c("00","01","10","11")
diag(Q_pos) <- -rowSums(Q_pos)

# neg Q matrix
x1=x2=n1=n2=associate_factor*sub_rate
y1=y2=m1=m2=sub_rate
Q_neg = matrix(c(0,x1,x2,0,
                y1,0,0,y2,
                m1,0,0,m2,
                0,n1,n2,0),
                4,4,byrow=TRUE)
rownames(Q_neg)<-colnames(Q_neg)<-c("00","01","10","11")
diag(Q_neg) <- -rowSums(Q_neg)

# simulate phenotype from 4Q matrix
tt<-sim.history(tree,t(Q_pos),message=FALSE)
tt1<-mergeMappedStates(tt,c("00","01"),"0")
tt1<-mergeMappedStates(tt1,c("10","11"),"1")
tt2<-mergeMappedStates(tt,c("00","10"),"0")
tt2<-mergeMappedStates(tt2,c("01","11"),"1")

x<-getStates(tt1,"tips")
y<-getStates(tt2,"tips") # pheno tip states

# get pheno at tips
pheno=as.numeric(y=="1")
names(pheno) = names(y)

cat("Phenotype 1 percentage")
sum(pheno)/length(pheno)

# infer a phenotype Q matrix based on correlated Q matrix
pheno_Q<-matrix(c(0,x1+m2,y1+n2,0),2,2,byrow=TRUE)
rownames(pheno_Q)<-colnames(pheno_Q)<-c(0,1)
diag(pheno_Q) <- -rowSums(pheno_Q)

# # infer a genotype Q matrix based on correlated Q matrix
# geno_Q<-matrix(c(0,x2+y2,m1+n1,0),2,2,byrow=TRUE)
# rownames(geno_Q)<-colnames(geno_Q)<-c(0,1)
# diag(geno_Q) <- -rowSums(geno_Q)
# geno_Q

# simulate to get phenotype ancestral state probability
cl<-makeSOCKcluster(rep("localhost", sim_ncpus))
sim_res<-snow::clusterApply(cl,
    x=replicate(sim_ncpus, pheno, simplify=FALSE),
    fun=make.simmap,tree=tree,Q=pheno_Q,nsim=sim_times/sim_ncpus)

sim_res<-do.call(c, sim_res)
if(!("multiSimmap"%in%class(sim_res)))
    class(sim_res)<-c("multiSimmap",
        class(sim_res))
snow::stopCluster(cl)

# get each internal node's state in each tree
clus <- makeSOCKcluster(rep("localhost", all_ncpus))
node.states.list=snow::parLapply(cl=clus, x=sim_res, fun=getStates, type="nodes")
snow::stopCluster(clus)

# infer internal pheno prob
node.1.prob = rep(0, tree$Nnode) # the prob for pheno==1
for (i in 1:length(sim_res)) {
    for(j in 1:tree$Nnode){
        if(node.states.list[[i]][j]==1){
            node.1.prob[j]=node.1.prob[j]+1
        }
    }
}
node.1.prob = node.1.prob/length(sim_res)
names(node.1.prob) = names(node.states.list[[1]])

# tip pheno prob
tip.1.prob = rep(0, length(tree$tip.label)) # the prob for pheno==1
for (i in 1:length(tree$tip.label)) {
    tip.1.prob = pheno
}
names(tip.1.prob) = 1:length(tree$tip.label)

# combine prob
pheno.1.prob = c(tip.1.prob, node.1.prob)

### simulate associated genotypes ###
get_associated_geno <- function(tree, pheno.1.prob, Q){
    tip_ct = length(tree$tip.label)
    node_ct = length(tree$node.label)
    root = tip_ct+1 # root index

    geno_states = pheno_states = rep(NA, tip_ct+node_ct)
    names(geno_states)=names(pheno_states)=names(pheno.1.prob)

    geno_states[root] = sample(c(0,1), 1, prob=c(0.5,0.5)) # root geno
    pheno_states[root] = sample(c(0,1), 1, prob=c(1-pheno.1.prob[tip_ct+1],pheno.1.prob[tip_ct+1])) # root pheno

    for(e in 1:length(tree$edge.length)){
        from_node = tree$edge[e,1]
        to_node = tree$edge[e,2]
        from_index = which(names(geno_states)==from_node)
        to_index = which(names(geno_states)==to_node)

        from_state = paste0(geno_states[from_index], pheno_states[from_index])

        P = matexpo(Q*tree$edge.length[e])
        colnames(P)=rownames(P)=colnames(Q)

        pheno_1_prob = pheno.1.prob[which(names(pheno.1.prob)==to_node)]
        prob = P[which(rownames(P)==from_state),]
        prob_mod = modify_prob(prob, c(1-pheno_1_prob, pheno_1_prob))
        to_state = sample(colnames(P), 1, prob=prob_mod)
        geno_states[to_index] = as.numeric(substr(to_state, 1, 1))
        pheno_states[to_index] = as.numeric(substr(to_state, 2, 2))
    }
    geno_tip_states = geno_states[1:tip_ct]
    names(geno_tip_states) = tree$tip.label
    return(geno_tip_states)
}

pos_geno = list()
for(k in 1:N_pos_geno){
    geno_tip_states = get_associated_geno(tree, pheno.1.prob, Q_pos)
    pos_geno[[paste0("pos", k)]] = geno_tip_states
}

neg_geno = list()
for(k in 1:N_neg_geno){
    geno_tip_states = get_associated_geno(tree, pheno.1.prob, Q_neg)
    neg_geno[[paste0("neg", k)]] = geno_tip_states
}

pos_geno_df = do.call(rbind, pos_geno)
neg_geno_df = do.call(rbind, neg_geno)
asso_geno_df = rbind.data.frame(pos_geno_df, neg_geno_df)

### simulate non-associated genotypes ###
# coverage is sampled from a beta distribution
N_nonasso = N_all_geno - N_pos_geno - N_neg_geno
set.seed(7)
coverages = abs(jitter(rbetabinom.ab(N_nonasso, size=Ntips, shape1=0.007, shape2=0.75)))/Ntips
coverages[coverages >= 1]=1-0.0001
coverages[coverages < 0.0001]=0.0001

length(coverages)
summary(coverages)

non_asso_geno <- foreach(i=1:length(coverages), .combine='rbind') %dopar% {
    coverage = coverages[i]
    alpha = coverage*2*(sub_rate+associate_factor*sub_rate) # make it reach equilibrium faster
    beta = (1-coverage)*2*(sub_rate+associate_factor*sub_rate)

    geno_Q <-matrix(c(0,alpha,beta,0),2,2,byrow=TRUE)
    rownames(geno_Q)<-colnames(geno_Q)<-c(0,1)
    diag(geno_Q) <- -rowSums(geno_Q)
    # keep sampling until at least one species has this gene
    while(TRUE){
        g_t <- sim.history(tree, t(geno_Q), message=FALSE)
        if(sum(g_t$states==1)!=0){
            non_asso_tip_state = as.numeric(g_t$states)
            break
        } 
    }
    non_asso_tip_state
}
rownames(non_asso_geno) = paste0("g", 1:nrow(non_asso_geno))

# dump tree, trait/pheno and geno data
geno_data = rbind.data.frame(non_asso_geno, pos_geno_df, neg_geno_df) %>% rownames_to_column(var="orthoID")
write_tsv(geno_data, paste0(outdir, "/gene.tsv"))
pheno_data = data.frame(Tip=names(pheno), Status=pheno)
write_tsv(pheno_data, paste0(outdir, "/trait.tsv"))
write.tree(tree, paste0(outdir, "/tree.nwk"))
