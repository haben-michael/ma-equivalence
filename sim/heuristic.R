source('utils.R')
set.seed(1)
ns <- 10:20
tau2.null <- tau2.star <- .25
re.gamma.centered <- function(n, shape = 10) {
  rgamma(n, shape = shape, rate = sqrt(shape)) - sqrt(shape)
}
start <- Sys.time()
by.n <- sapply(ns, function(n) {
    cat('.')
    m <- round(n^1)
    stats <- replicate(n^2, {
        yvb <- dgp.logodds(n=n,m=m,tau2.star=tau2.star,rho=1.5,re.distr=re.gamma.centered)
        y <- yvb$y; v <- yvb$v; b <- yvb$b; sigma2 <- v/m
        test.stat.LR <- ma.LR(y,sigma2)['test.stat']
        test.stat.REML <- with(metafor::rma(yi=y,vi=sigma2,method='REML',control=list(stepadj=.5)), zval)
        c(LR=test.stat.LR, REML=test.stat.REML)
    })
    ecdf.LR <- ecdf(pchisq(stats['LR',],df=1,lower=FALSE))
    ecdf.REML <- ecdf(2*pnorm(abs(stats['REML',]),lower=FALSE))
    c(LR=max(abs(ecdf.LR(knots(ecdf.LR)) - knots(ecdf.LR))),
    REML=max(abs(ecdf.REML(knots(ecdf.REML)) - knots(ecdf.REML))))   
})
by.n <- simplify2array(by.n)
print(Sys.time()-start)
