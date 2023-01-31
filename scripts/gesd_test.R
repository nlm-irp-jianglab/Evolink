library(PMCMRplus) # need to conda install r-gmp and r-mpfr

grubbsTest <- function (x, m=0, std=NULL, alternative = c("two.sided", "greater", "less"))
{
    if (!is.numeric(x)) {
        stop("'x' must be a numeric vector.")
    }
    alternative <- match.arg(alternative)
    data.name <- deparse(substitute(x))
    x <- na.omit(x)
    Mean <- m
    S <- std
    n <- length(x)
    if (alternative == "two.sided") {
        MxAbs <- max(abs(x - Mean))
        G <- MxAbs/S
        pval <- min(1, 2 * pgrubbs(G, n, lower.tail = FALSE))
        i <- which(MxAbs == abs(x - Mean))
        val <- x[i]
    }
    else if (alternative == "greater") {
        val <- max(x)
        G <- (val - Mean)/S
        pval <- pgrubbs(G, n, lower.tail = FALSE)
        i <- which(val == x)
    }
    else {
        val <- min(x)
        G <- (Mean - val)/S
        pval <- pgrubbs(G, n, lower.tail = FALSE)
        i <- which(val == x)
    }
    names(val) <- NULL
    names(i) <- NULL
    ans <- list(method = "Grubbs single outlier test", alternative = alternative,
        statistic = c(G = G), parameter = c(df = n - 2), p.value = pval,
        estimate = c(c(i = i), c(value = val)), data.name = data.name)
    class(ans) <- "htest"
    return(ans)
}

gesdTest <- function(x, alt="two.sided"){
    x <- na.omit(x)
    n <- length(x)
    oldx <- x
    ix <- rep(NA, n)
    PVAL <- rep(NA, n)
    R <- rep(NA, n)
    std = sd(oldx)
    
    ## repeated single outlier Grubb's test
    while(TRUE){
        if(length(x)==0) break
        out <- grubbsTest(x, m=0, std=std, alternative = alt)
        cur_indices = out$estimate[grep("i", names(out$estimate))] # return the index in the current vector
        values = out$estimate[grep("value", names(out$estimate))]
        indices = which(oldx %in% values)
        ix[indices] <- indices
        PVAL[indices] <- out$p.value
        R[indices] <- out$statistic
        x <- x[-cur_indices]
    }

    ans <- list(method = "GESD multiple outlier test",
                statistic = R,
                p.value = PVAL,
                ix = ix,
                alternative = "two.sided")
    class(ans) <- "gesdTest"
    return(ans)
}