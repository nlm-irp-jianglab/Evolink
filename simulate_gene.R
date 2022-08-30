suppressPackageStartupMessages({
    suppressWarnings(library(phangorn))
    suppressWarnings(library(geiger))
})

.is.integer0 <- function(x){
  #### function from treeWAS (https://github.com/caitiecollins/treeWAS) ####
  is.integer(x) && length(x) == 0L
}

get.fitch.n.mts <- function(x, tree, snps=NULL){
  #### function from treeWAS (https://github.com/caitiecollins/treeWAS) ####
  ## Re-coding snps as x (to allow for phen/vectors).
  ## --> snps now deprecated:
  X <- NULL
  if(!missing(x)){
    X <- x
    if(!is.null(snps) & !is.null(x)){
      warning("As 'x' is specified, we ignore the 'snps' argument. \n
              (In get.fitch.n.mts the 'snps' argument has now been replaced by an argument named 'x'.)")
    }
  }else{
    if(!is.null(snps)){
      X <- snps
    }
  }
  ## If ONE of x or snps was specified, continue; else, stop:
  if(!is.null(X)){
    x <- X
  }else{
    stop("'x' must be specified.")
  }

  ## checks
  # if(!is.matrix(x)) stop("x must be a matrix.") ## no it mustn't... works for phen too.
  # levs <- unique(as.vector(unlist(x)))
  levs <- unique(as.vector(as.matrix(x)))
  # if(any(is.na(levs))) levs <- levs[-which(is.na(levs))]
  if((!is.numeric(x) & !is.logical(x)) | length(levs[!is.na(levs)])!=2){
    stop("x must be a numeric matrix or vector, with two unique values, excluding NAs
         (though we recommend that NAs be in the minority for each column).")
  }
  if(any(is.na(levs))){
    if(is.matrix(x)){
      nnas <- sapply(c(1:ncol(x)), function(e) length(which(is.na(x[,e])))/nrow(x))
      toRemove <- which(nnas > 0.5)
      if(length(toRemove) > 0){
        cat(length(toRemove), "snps columns are over 50% NAs.
            You may want to consider removing these columns as they are unlikely to be significant
            and can generate inappropriate inferences during ancestral state reconstruction.")
      }
    }else{
      nnas <- length(which(is.na(x)))/length(x)
      # toRemove <- which(nnas > 0.5)
      if(nnas > 0.5){
        cat("x is over 50% NAs.
            This may generate inappropriate inferences during ancestral state reconstruction.")
      }
    }
  }

  x.levels <- sort(levs, na.last = TRUE)
  ## returns only unique patterns...
  x.phyDat <- as.phyDat(as.matrix(x),
                           type="USER", levels=x.levels)
  ## get index of all original x columns to map to unique pattern
  index <- attr(x.phyDat, "index")
  
  ## get parsimony score for all unique patterns in x
  ## NB: For fitch.phangorn, x data must be of class phyDat
  fitch.unique <- phangorn::fitch(tree, x.phyDat, site="site")
  # table(fitch.unique)

  ## get score for all original sites
  fitch.complete <- fitch.unique[index]
  return(fitch.complete)
}

#### get background by simulation ####
gene.sim <- function(tree = NULL, # must be given and must be rooted: if(!is.rooted(tree)) tree <- midpoint(tree)
                     gene = NULL, # must be given, row=genes, col=species
                     sim.fold = 10, # how many times of number the genotypes should be simulated 
                     seed = 1){

  # @tree: must be given and must be rooted: if(!is.rooted(tree)) tree <- midpoint(tree)
  # @gene: must be given, row=genes, col=species
  # @sim.fold: default=10, how many times we will simulate on the observed genotypes
  # @seed: set seed to make results reproducable

  # suppress warnings
  defaultW <- getOption("warn")
  options(warn = -1)

  # cls = makeCluster(clstr, type = "SOCK")

  gene = t(gene) # gene matrix row must be tree tips
  n.subs <- get.fitch.n.mts(x=gene, tree=tree)
  n.sim.gene = ncol(gene)*sim.fold #number of simulated genotypes

  ## Always work with trees in "pruningwise" order:
  tree <- reorder.phylo(tree, order="pruningwise")

  ## Get sequence from lowest ("root", Nterm+1) to highest ancestral node:
  ix <- c(min(tree$edge[,1]):max(tree$edge[,1]))
  ## Get for loop index of rows in tree$edge[,1], in pairs (if binary tree), from lowest to highest:
  x <- as.vector(unlist(sapply(c(1:length(ix)), function(e) which(tree$edge[,1] == ix[e]))))

  ## GET MUTATIONS' branch & loci ##
  n.ind <- min(tree$edge[,1])-1 # tree$Nnode+1

  gen.size <- n.sim.gene
  edges <- tree$edge

  ## Simulate genotype for root individual: ##
  if(!is.null(seed)) set.seed(seed)
  gen.root <- sample(c(TRUE, FALSE), gen.size, replace=TRUE) # assign binary state for each gen.size

  ## get the sum of all branch lengths in the tree:
  time.total <- sum(tree$edge.length)

  ## make dummy variables in which to store the resulting n.mts variables:
  n.mts <- new.nt <- NA

  #### handle nsubs ####
  ## a vector (containing a distribution)
  ## --> use this distribution to define n.subs-per-site
  ## if a distribution is provided by the user,
  ## we use this to determine the number of substitutions
  ## to occur at what proportion of sites (note that
  ## we may not be simulating the same number of sites)

  dist <- n.subs

  ## get dist.prop, a distribution containing the counts
  ## of the genes to be simulated that will have
  ## i many substitutions
  dist.sum <- sum(dist) # should be the number of genes
  dist.prop <- round((dist/dist.sum)*gen.size) # gen.size=n.genes.sim, but dist.sum=n.genes, round() to have integers

  ## check that these counts sum to gen.size,
  ## else add the remainder to the largest n.subs count
  ## (OR should we just add these to the n.subs=1 set ??? ###
  ## likely to be the same thing, but could not be...)
  # after round(), if have more/less, then reduce/increase the top max element values to make the sum(dist.prop) == gen.size
  dist.prop.bak = dist.prop
  
  while(sum(dist.prop.bak) != gen.size){
      m <- which.max(dist.prop)
      if(sum(dist.prop.bak) > gen.size){delta = 1}else{delta = -1}
      if((dist.prop.bak[m] - delta) >= 0){
        dist.prop.bak[m] = dist.prop.bak[m] - delta
        dist.prop = dist.prop[-m]
      }
	    if(m==1){
        dist.prop = dist.prop.bak
      }
  }

  dist.prop = dist.prop.bak

  ## get rid of useless trailing 0s, if the last element is 0, cut it off from dist to shorten the dist
  while(dist.prop[length(dist.prop)] == 0){
    dist.prop <- dist.prop[c(1:(length(dist.prop)-1))]
  }

  ## make n.mts, a vector of length ncol(genes)
  n.mts <- rep(1111, gen.size)
  loci.available <- c(1:gen.size)
  ## assign dist.prop[i] elements of n.mts
  ## to be the same as the n.subs
  ## indicated by i, the given element of dist.prop: e.g. 
  for(j in 1:length(dist.prop)){
    ## provided there are not 0 sites to have this number of substitutions...
    if(dist.prop[j] > 0){
      if(length(loci.available) > 1){
        if(!is.null(seed)) set.seed(seed)
        ## assign dist.prop[i] elements of n.mts to be i
        loci.selected <- sample(loci.available, dist.prop[j], replace = FALSE) # choose dist.prop[j] sites from loci.available
        loci.available <- loci.available[-which(loci.available %in% loci.selected)] # update loci.available
      }else{
        ## if there is only 1 (the last) loci available,
        ## we select this one:
        loci.selected <- loci.available
      }
      n.mts[loci.selected] <- j
    }
  }
  # n.mts in the end is a vector with gen.size/n.sim.genes sites whose n.subs are sampled from the dist.prop.
  rm(dist)
  rm(dist.prop)
  rm(dist.sum)

  ## get non-associated genes ##

  ## for each site, draw the branches to which
  ## you will assign the mts for this site
  ## (~ branch length):

  l.edge <- length(tree$edge.length)
  # default FALSE for all edges
  null.vect <- rep(FALSE, l.edge)
  
  if(max(n.mts) > l.edge){
    if(!is.null(seed)) set.seed(seed)
    genes.loci <- NULL
      for(e in 1:length(n.mts)){
        genes.loci <- list()
        repTF <- FALSE
        if(n.mts[e] > l.edge) repTF <- TRUE
        genes.loci[[e]] <- replace(null.vect, sample(c(1:l.edge),
                                                    n.mts[e],
                                                    replace=repTF,
                                                    prob=tree$edge.length), TRUE)

        genes.loci <- t(do.call(rbind, genes.loci))
      }

  }else{ # here
    if(!is.null(seed)) set.seed(seed)
    # remember n.mts is a vector for all n.sim.genes sites' mutation count
    # for each site, we randomly assign n.mts[e] mutations to edges acrding to the tree length (namely longer edge is more likely to have a mutation)
    # genes.loci is a M*N matrix, with M=edge number, N=n.sim.genes
    genes.loci <- sapply(c(1:length(n.mts)),
                        function(e)
                          replace(null.vect,
                          sample(c(1:l.edge),
                                 n.mts[e],
                                 replace=FALSE,
                                 prob=tree$edge.length), TRUE))
    # clusterExport(cls, c("null.vect", "l.edge", "tree", "n.mts"))
    # mut_edge = function(e){
    #     x=replace(null.vect,
    #             sample(c(1:l.edge),
    #             n.mts[e],
    #             replace=FALSE,
    #             prob=tree$edge.length), 
    #             TRUE)
    #     return(x)
    # }
    # genes.loci <- parSapply(cls, c(1:length(n.mts)), mut_edge)
  }

  ## rearrange genes.loci s.t it becomes a
  ## list of length tree$edge.length,
  ## each element of which contains the
  ## **n.genes.sim** locations of the mutations that will
  ## occur on that branch
  genes.loci <- sapply(c(1:nrow(genes.loci)),
                       function(e) which(genes.loci[e,] == TRUE))                       

  ## get the node names for all individuals (terminal and internal) ==> all node names (tips and internal nodes)
  all.inds <- sort(unique(as.vector(unlist(tree$edge))))
  # we will store the output in a list called genes:
  genes <- list()
  ## we start w all inds having same genotype as root:
  # genes is **a list of tree nodes (e.g. 198) long with each element a vector of genes' binary states (e.g. 200)** (initially assigned with root states)**
  genes[all.inds] <- rep(list(gen.root), length(all.inds))

  ## store replacement nts in list new.nts:
  new.nts <- list()
  ## distinguish btw list of loci and unique list
  genes.loci.ori <- genes.loci
  ## will need to treat repeat loci differently
  genes.loci.unique <- lapply(genes.loci, unique)

  ## For Loop to get new nts ##
  for(i in x){
    # iterate edge index

    ## for all genes other than root, we mutate the
    ## genome of the node preceding it, according to genes.loci.
    ## Draw new nts for each locus selected for mutation:
    if(!.is.integer0(genes.loci.unique[[i]])){
      # tree$edge[i,1] is "from_node" and genes.loci.unique[[i]] indicates the genes.loci should mutate now and "!" is used to mutate
      new.nts[[i]] <- !genes[[tree$edge[i,1]]][genes.loci.unique[[i]]] # for this edge, new.nts is updated

      ## if any loci are selected for multiple mutations (namely genes.loci.ori[[i]] has redundant loci.)
      ## within their given branch length:
      if(length(genes.loci.ori[[i]]) != length(genes.loci.unique[[i]])){
        ## identify which loci are repeaters; find the redundant loci
        repeats <- table(genes.loci.ori[[i]])[which(table(genes.loci.ori[[i]])!=1)] 
        ## how many times they repeat
        n.reps <- repeats - 1
        ## the positions of these loci in the vector of genes loci
        toRepeat <- which(genes.loci.unique[[i]] %in% names(repeats))
        ## run chain of re-sampling to end in our new nt for repeater loci:
        foo <- list()
        for(j in 1:length(toRepeat)){
          foo[[j]] <- new.nts[[i]][toRepeat[j]]
          for(k in 1:n.reps[j]){
            if(k==1){
              foo[[j]][k] <- !foo[[j]][1]

            }else{
              foo[[j]][k] <- !foo[[j]][k-1]
            }
          }
          ## retain only the last nt selected
          out <- sapply(c(1:length(foo)),
                        function(e) foo[[e]][length(foo[[e]])])
        }
        ## for the loci with repeated mts, replace these positions
        ## in new.nts with the corresponding elements of out, above.
        new.nts[[i]][toRepeat] <- out
      } # end of if statement for repeaters

      ## update ancestral genotype with new.nts:
      temp <- genes[[tree$edge[i,1]]] # get all the states for gen.size loci at this from_node
      temp[genes.loci.unique[[i]]] <- new.nts[[i]] # get loci states if they have been updated
      genes[[tree$edge[i,2]]] <- temp # pass the new states to the to_node in genes

    }else{
      ## if no mts occur on branch, set genotype of
      ## downstream individual to be equal to ancestor's
      genes[[tree$edge[i,2]]] <- genes[[tree$edge[i,1]]] # if no change, just copy the states from from_node to to_node
    }
  } # end of for loop selecting new nts at mutator loci
  # Eventually, genes records each node's sim.genes states after mutation.

  ## CHECK IF ALL LOCI ARE POLYMORPHIC (|polyThres) ##

  ## temporarily assemble non-associated loci into matrix:
  # temp.ori <- do.call("rbind", genes)
  temp <- do.call("rbind", genes) # convert genes from a list to N_node*n.sim.genes matrix

  ## keep only rows containing terminal individuals:
  # temp.ori <- temp.ori[1:n.ind, ]
  temp <- temp[1:n.ind, ] # only get the tip states => n.tip*n.sim.genes matrix


  ## identify n.minor.allele required to meet polyThres:
  polyThres <- 0.001
  n.min <- n.ind*polyThres

  ## make a list of any NON-polymorphic loci:
  csum <- colSums(temp) # the positive states for each loci/genotype
  toRepeat <- which(csum < n.min | csum > (nrow(temp) - n.min)) # n.min=minimal number of 1s; "nrow(temp)-n.min"=minimal number of 0s

  ## REPLACE ANY NON-POLYMORPHIC LOCI & GENERATE ASSOCIATED genes ##

  ## NEW while loop ##

  counter <- 0
  if(!is.null(seed)) seed.i <- seed

  while(length(toRepeat) > 0){
    # print("LENGTH toREPEAT"); print(length(toRepeat))
    ## for each site, draw the branches to which
    ## you will assign the mts for this site
    ## (~ branch length):

    ## Get vector of FALSEs of length tree$edge.length:
    null.vect <- rep(FALSE, length(tree$edge.length))

    if(!is.null(seed)){
      seed.i <- seed.i+1
      set.seed(seed.i)
    }

    if(max(n.mts[toRepeat]) > length(tree$edge.length)){rep_tf = TRUE}else{rep_tf = FALSE}
    genes.loci <- sapply(c(1:length(n.mts[toRepeat])),
                        function(e)
                          replace(null.vect,
                                  sample(c(1:length(tree$edge.length)),
                                          n.mts[toRepeat][e],
                                          replace=rep_tf,
                                          prob=tree$edge.length), TRUE))
    # clusterExport(cls, c("null.vect", "rep_tf", "toRepeat", "n.mts", "tree"))
    # mut_edge = function(e){
    #     x=replace(null.vect,
    #             sample(c(1:length(tree$edge.length)),
    #                     n.mts[toRepeat][e],
    #                     replace=rep_tf,
    #                     prob=tree$edge.length), 
    #             TRUE)
    #     return(x)
    # }
    # genes.loci <- parSapply(cls, c(1:length(n.mts[toRepeat])), mut_edge)

    ## rearrange genes.loci s.t it becomes a
    ## list of length tree$edge.length,
    ## each element of which contains the
    ## locations of the mutations that will
    ## occur on that branch
    genes.loci <- sapply(c(1:nrow(genes.loci)),
                        function(e) which(genes.loci[e,] == TRUE))


    ## get the node names for all individuals (terminal and internal)
    all.inds <- sort(unique(as.vector(unlist(tree$edge))))
    for(i in 1:length(all.inds)){
      genes[[all.inds[i]]][toRepeat] <- gen.root[toRepeat]
    }

    ## store replacement nts in list new.nts:
    new.nts <- list()
    ## distinguish btw list of loci and unique list
    genes.loci.ori <- genes.loci
    ## will need to treat repeat loci differently...
    genes.loci.unique <- lapply(genes.loci, unique) # (identical to genes.loci?)

    ## for loop to get new nts ##

    for(i in x){
      ## for all genes other than root, we mutate the
      ## genome of the node preceding it, according to genes.loci.
      ## Draw new nts for each locus selected for mutation:
      if(!.is.integer0(genes.loci.unique[[i]])){
        new.nts[[i]] <- !genes[[tree$edge[i,1]]][toRepeat][genes.loci.unique[[i]]]

        ## if any loci are selected for multiple mutations
        ## within their given branch length:
        if(length(genes.loci.ori[[i]]) != length(genes.loci.unique[[i]])){
          ## identify which loci are repeaters
          repeats <- table(genes.loci.ori[[i]])[which(table(genes.loci.ori[[i]])!=1)]
          ## how many times they repeat
          n.reps <- repeats - 1
          ## the positions of these loci in the vector of genes loci
          toRep <- which(genes.loci.unique[[i]] %in% names(repeats))
          ## run chain of re-sampling to end in our new nt for repeater loci:
          foo <- list()
          for(j in 1:length(toRep)){
            foo[[j]] <- new.nts[[i]][toRep[j]]
            for(k in 1:n.reps[j]){
              if(k==1){
                foo[[j]][k] <- !foo[[j]][1]
              }else{
                foo[[j]][k] <- !foo[[j]][k-1]
              }
            }
            ## retain only the last nt selected
            out <- sapply(c(1:length(foo)),
                          function(e) foo[[e]][length(foo[[e]])])
          }
          ## for the loci with repeated mts, replace these positions
          ## in new.nts with the corresponding elements of out, above.
          new.nts[[i]][toRep] <- out
        } # end of if statement for repeaters

        ## update ancestral genotype with new.nts:
        temp1 <- genes[[tree$edge[i,1]]][toRepeat]
        temp1[genes.loci.unique[[i]]] <- new.nts[[i]]
        genes[[tree$edge[i,2]]][toRepeat] <- temp1

      }else{
        ## if no mts occur on branch, set genotype of
        ## downstream individual to be equal to ancestor's
        genes[[tree$edge[i,2]]][toRepeat] <- genes[[tree$edge[i,1]]][toRepeat]
      }
    } # end of for loop selecting new nts at mutator loci


    ## temporarily assemble non-associated loci into matrix:
    # temp.new <- do.call("rbind", genes)
    temp <- do.call("rbind", genes)

    ## keep only rows containing tip nodes
    temp <- temp[1:n.ind, ]

    ## CHECK IF ALL LOCI ARE POLYMORPHIC (|polyThres)

    ## identify n.minor.allele required to meet polyThres:
    polyThres <- 0.01
    n.min <- n.ind*polyThres

    ## make a list of any NON-polymorphic loci:
    toRepeat.ori <- toRepeat
    temp.toRepeat <- temp[, toRepeat.ori]

    ## make a vector of any NON-polymorphic loci:
    if(!is.matrix(temp.toRepeat)){
      csum <- sum(temp.toRepeat)
      if(csum < n.min | csum > (length(temp.toRepeat) - n.min)){
        toRepeat <- toRepeat.ori
      }else{
        toRepeat <- NULL
      }
    }else{
      ## if temp.toRepeat is a true matrix:
      if(ncol(temp.toRepeat) > 0){
        csum <- colSums(temp.toRepeat)
        toRepeat <- toRepeat.ori[which(csum < n.min | csum > (nrow(temp.toRepeat) - n.min))]
      }else{
        csum <- sum(temp.toRepeat)
        if(csum < n.min | csum > (length(temp.toRepeat) - n.min)){
          toRepeat <- toRepeat.ori
        }else{
          toRepeat <- NULL
        }
      }
    }

    counter <- counter+1
    # print("COUNTER"); print(counter)
    # print("toRepeat"); print(length(toRepeat))
  }

  # Prepare output
  colnames(temp) <- c(1:ncol(temp))

  rm(genes)
  sim_geno <- temp
  rm(temp)

  # stopCluster(cls)

  ## keep only rows containing terminal individuals:
  sim_geno <- sim_geno[1:n.ind, ]

  ## assign rownames to match tree$tip.label:
  if(!is.null(tree$tip.label)){
    if(length(tree$tip.label) == nrow(sim_geno)){
      rownames(sim_geno) <- tree$tip.label
    }else{
      rownames(sim_geno) <- 1:nrow(sim_geno)
      warning("The length of tree$tip.label was not equal to nrow(genes) being simulated;
            rownames(genes) have been set to 1:N and will not match tree$tip.label.")
    }
  }

  ## generate column names:
  colnames(sim_geno) <- paste0("s", 1:ncol(sim_geno))

  ## convert bool to numeric without changing rownames/colnames
  sim_geno[] = sapply(sim_geno, as.integer)

  ## get results: ##
  return(data.frame(t(sim_geno)))
}
