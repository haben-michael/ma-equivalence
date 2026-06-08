args <- commandArgs(trailingOnly=TRUE)
m.exp <- as.numeric(args[1])

source('utils.R')
ns <- 10:25
corr.yv <- .4
bias.y <- 1
tau2.null <- tau2.star <- 1
estimator.strings <- c(DL='DL',HE='HE',HS='HS',SJ='SJ',REML='REML',PM='PM')
start <- Sys.time()
by.n <- lapply(ns, function(n) {
    m <- n^m.exp
    stats <- replicate(n^3, expr={
        yv <- dgp.gamma(n,var.y=tau2.star,var.v=1,corr.yv=corr.yv)
        y <- yv$y+bias.y/m; v <- yv$v; sigma2 <- v/m
        mle.full <- mle.ht(y=y,sigma2=sigma2)
        mu.mle.full <- mle.full['mu']; tau2.mle.full <- mle.full['tau2']
        IV.tau2s <- sapply(estimator.strings,function(str)
            metafor::rma(yi=y,vi=sigma2,method=str,control=list(maxiter=1000,stepadj=.5))$tau2)
        tau2.hats <- c(IV.tau2s,naive=(n-1)/n*var(y))
        tau2.hats - tau2.mle.full
    })
    apply(stats,1,function(x)median(abs(x)))
})
by.n <- simplify2array(by.n)


filename <- paste0('save',as.integer(abs(rnorm(1))*1e8),'.RData')
save(m.exp,corr.yv,by.n=by.n,ns=ns,file=filename)


