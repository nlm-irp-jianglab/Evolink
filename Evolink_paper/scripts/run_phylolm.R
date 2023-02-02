library(phylolm)
library(tidyverse)
library(optparse)
library(foreach)
library(doParallel)
registerDoParallel(56)
set.seed(1234567)
 
option_list = list(
  make_option(c("-n", "--tree"), type="character", default=NULL, 
              help="ultrametric tree file path", metavar="character"),
  make_option(c("-g", "--patable"), type="character", default=NULL, 
              help="gene PA matrix file path", metavar="character"),
  make_option(c("-t", "--trait"), type="character", default=NULL, 
              help="feature/trait", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="test.out", 
              help="output path", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read phylum time tree
tree_path = opt$tree
tre = read.tree(tree_path) 
taxa = sort(tre$tip.label) # sort taxa

# presence/absence gene matrix
patable = read_tsv(opt$patable, show_col_types = FALSE)
PA_df = patable %>% column_to_rownames(var="gene") %>% t() %>% as.data.frame() %>% rownames_to_column(var="genomeID")
x = PA_df %>% filter(genomeID %in% taxa) %>% arrange(genomeID)

# feature/trait for genomes
habitat_path = opt$trait
habitat_df = read_tsv(habitat_path)
y <- habitat_df %>% filter(Tip %in% taxa) %>% arrange(Tip)

get_btol <- function(my_trait01, my_predictor, nonNA_taxa){
    test_dat = data.frame(trait01 = my_trait01, predictor = my_predictor, row.names=c(nonNA_taxa))
    mf = model.frame(trait01~predictor, data=test_dat)
    mf = mf[which(rownames(mf) %in% nonNA_taxa),,drop=F]
    X = model.matrix(attr(mf, "terms"), data=mf)
    fit = glm(my_trait01~X-1, family=binomial)
    startB = fit$coefficients
    mybtol <- max(ceiling(max(abs(X%*%startB))), 10)
    return(mybtol)
}

# fit phyloglm regression model
start=2
end <- length(colnames(x))

res <- foreach(i=start:end, .combine=rbind) %dopar% {
    orthoID = colnames(PA_df)[i]
    # print(paste0(i, ":", orthoID))

    # customize btol for each gene
    nonNA_taxa <- y %>% filter(!is.na(Status)) %>% .[["Tip"]]
    my_trait01 <- y %>% filter(Tip %in% nonNA_taxa) %>% .[["Status"]]
    my_predictor <- x %>% filter(genomeID %in% nonNA_taxa) %>% .[[orthoID]]
    mybtol <- get_btol(my_trait01, my_predictor, nonNA_taxa)

    if(var(my_predictor)==0){
        pred_par <- c(NA, NA)
        names(pred_par) <- c("Estimate", "p.value")
        int_par <- c(NA, NA)
        names(int_par) <- c("Estimate", "p.value")
        aic = NA
        logLik = NA
        penlogLik = NA
    }else{
        dat = data.frame(trait01 = y$Status, predictor = x[[orthoID]], row.names=c(taxa))
        fit <- tryCatch(
            phyloglm(trait01~predictor, data=dat, phy=tre, btol=mybtol),
            error = function(e) {
                warning(paste(e))
                c(NA)
            })
        if(is.na(fit)){
            # if error occurs, we set btol down to 30 to make it work
            fit <- phyloglm(trait01~predictor, data=dat, phy=tre, btol=30)
            # pred_par <- c(NULL, NULL)
            # names(pred_par) <- c("Estimate", "p.value")
            # int_par <- c(NULL, NULL)
            # names(int_par) <- c("Estimate", "p.value")
            # aic = NULL
            # logLik = NULL
            # penlogLik = NULL
        }
        pred_par = summary(fit)$coefficients["predictor", c("Estimate", "p.value")]
        int_par = summary(fit)$coefficients["(Intercept)", c("Estimate", "p.value")]
        aic = summary(fit)$aic
        logLik = summary(fit)$logLik
        penlogLik = summary(fit)$penlogLik
    }
    c(orthoID, pred_par, int_par, aic, logLik, penlogLik)
}

res <- res %>% as.data.frame()
colnames(res) <- c("orthoID", "pred.Estimate", "pred.Pvalue", 
                   "int.Estimate", "int.Pvalue", 
                   "AIC", "logLik", "penlogLik")

# write to file
write_tsv(res, opt$output)
