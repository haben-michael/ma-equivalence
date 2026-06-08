args <- commandArgs(trailingOnly=TRUE)
m.exp <- as.numeric(args[1])
p.C.var <- .1^2

source('utils.R')
ns <- 10:20
tau2.null <- tau2.star <- .25
estimator.strings <- c(DL='DL',HE='HE',HS='HS',SJ='SJ',REML='REML',PM='PM')
re.gamma.centered <- function(n, shape = 10) {
  rgamma(n, shape = shape, rate = sqrt(shape)) - sqrt(shape)
}
start <- Sys.time()
by.n <- sapply(ns, function(n) {
    cat('.')
    m <- round(n^m.exp)
    stats <- replicate(n^3, {
        yvb <- dgp.logodds(n=n,m=m,tau2.star=tau2.star,rho=1.5,re.distr=re.gamma.centered)
        y <- yvb$y; v <- yvb$v; b <- yvb$b; sigma2 <- v/m
        mle.full <- mle.ht(y=y,sigma2=sigma2)
        mu.mle.full <- mle.full['mu']; tau2.mle.full <- mle.full['tau2']
        test.stats.IV <- sapply(estimator.strings,function(str)
            metafor::rma(yi=y,vi=sigma2,method=str,control=list(maxiter=1000,stepadj=.5))$zval)
        tau2.naive <- (n-1)/n*var(y)
        w <- 1/(tau2.naive+sigma2)
        test.stat.naive <- sum(y*w) / sqrt(sum(w))
        test.stats.IV <- c(test.stats.IV,naive=test.stat.naive)
        test.stat.LR <- ma.LR(y,sigma2)['test.stat']
        c(test.stats.IV,LR=test.stat.LR)
    })
    ecdf.LR <- ecdf(pchisq(stats['LR',],df=1,lower=FALSE))
    ecdfs.IV <- apply(subset(t(stats),select=-LR), 2, function(IV.stats)ecdf(2*pnorm(abs(IV.stats),lower=FALSE)),simplify=FALSE)
    ks.LR <- max(abs(ecdf.LR(knots(ecdf.LR)) - knots(ecdf.LR)))
    ks.IV <- sapply(ecdfs.IV, function(ecdf.IV) max(abs(ecdf.IV(knots(ecdf.IV)) - knots(ecdf.IV))))
    c(LR=ks.LR, ks.IV)
})
by.n <- simplify2array(by.n)

filename <- paste0('save',as.integer(abs(rnorm(1))*1e8),'.RData')
save(m.exp,p.C.mean,by.n=by.n,ns=ns,file=filename)


