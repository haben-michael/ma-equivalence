mle.ht <- function(mu.start=NULL,tau2.start=NULL,y,sigma2,B=15,
                   mu.update=NULL,tau2.update=NULL  ) {
    n <- length(y)
    if(is.null(mu.start)) {
        w.fe <- 1/sigma2
        mu.start <- y %*% w.fe / sum(w.fe)
    }
    if(is.null(tau2.start)) {
        w.fe <- 1/sigma2
        mu.hat.fe <- y %*% w.fe / sum(w.fe)
        Q <- w.fe%*%(y-c(mu.hat.fe))^2
        tau2.start <- max(0, (Q - (n-1)) / (sum(w.fe)-sum(w.fe^2)/sum(w.fe)))
    }
    if(is.null(mu.update)) mu.update <- function(tau2)sum(y / (sigma2+tau2)) / sum(1/(sigma2+tau2))
    if(is.null(tau2.update)) tau2.update <- function(mu,tau2)max(0,sum( ((y-mu)^2 - sigma2)/(sigma2+tau2)^2) / sum(1/(sigma2+tau2)^2))
    mu.mle <- mu.start; tau2.mle <- tau2.start
    for(i in 1:B) {
        ## mu.mle <- sum(y / (sigma2+tau2.mle)) / sum(1/(sigma2+tau2.mle))
        mu.mle <- mu.update(tau2.mle)
        ## tau2.mle <- sum( ((y-mu.mle)^2 - sigma2)/(sigma2+tau2.mle)^2) / sum(1/(sigma2+tau2.mle)^2)
        ## tau2.mle <- max(0,tau2.mle)
        tau2.mle <- tau2.update(mu.mle,tau2.mle)
    }
    return(c(mu=mu.mle,tau2=tau2.mle))
}


ma.LR <- function(y,sigma2,mu.null=0) {
    n <- length(y)    
    mle.full <- mle.ht(y=y,sigma2=sigma2)
    mu.mle.full <- mle.full['mu']; tau2.mle.full <- mle.full['tau2']
    mle.redu <- mle.ht(y=y,sigma2=sigma2,mu.update=function(tau2)mu.null)
    mu.mle.redu <- mle.redu['mu']; tau2.mle.redu <- mle.redu['tau2']
    L <- function(mu,tau2) -1/2*sum(log(2*pi*(sigma2+tau2))) - sum((y-mu)^2/2/(sigma2+tau2))
    test.stat <- 2*(L(mu.mle.full,tau2.mle.full) - L(mu.mle.redu,tau2.mle.redu))
    test.stat.bartlett <- test.stat / (1+3/2/n)
    return(structure(unname(c(test.stat, test.stat.bartlett)), names=c('test.stat','test.stat.bartlett')))
}



dgp.gamma <- function(n, var.y, var.v, corr.yv, rate = 1/sqrt(var.v)) {
  S1 <- var.y*(rate^2) 
  S2 <- var.v*(rate^2) 
  a3 <- corr.yv*sqrt(S1*S2)
  a1 <- S1 - a3
  a2 <- S2 - a3
  if (any(c(a1, a2, a3) <= 0)) stop()
  
  x1 <- rgamma(n, shape=a1, rate=rate)
  x2 <- rgamma(n, shape=a2, rate=rate)
  x3 <- rgamma(n, shape=a3, rate=rate)  
  y.unshifted <- x1 + x3
  y <- y.unshifted - (var.y*rate)
  v <- x2 + x3 
  return(list(y=y, v=v))
}

dgp.logodds <- function(n,m=n,tau2.star,mu.null=0,rho=1,re.distr=rnorm) { 
    b <- re.distr(n)*sqrt(tau2.star)+mu.null
    p.C <- plogis(b)
    m.C <- m; m.T <- round(rho*m)
    p.T <- plogis(qlogis(p.C)+b)
    x.C <- x.C.raw <- rbinom(n,m.C,p.C)
    x.T <- x.T.raw <- rbinom(n,m.T,p.T)
    a <- x.T.raw + 0.5      
    b <- m.T - x.T.raw + 0.5    
    c <- x.C.raw + 0.5          
    d <- m.C - x.C.raw + 0.5    
    y <- log((a / b) / (c / d))
    se2 <- 1/a + 1/b + 1/c + 1/d
    v <- mean(c(m.C, m.T)) * se2
    return(list(y=y,v=v,b=b))
}
