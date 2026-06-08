## plot routine for heuristic figure

plot(ns,abs(by.n['LR',]),type='l',ylim=range(abs(by.n)),xlab='n',ylab='distance to reference distribution')
lines(ns,abs(by.n['REML',]),lty=2)
legend('bottomleft',lty=1:2,legend=c('LR test','IV test'))
tikzDevice::tikz('./figs/heuristic.tex', width = 3.25, height = 2.75, pointsize = 8)
op <- par(mai=c(0.38, 0.48, 0.04, 0.03), mgp=c(1.5, 0.45, 0), cex.axis=0.60, cex.lab = 0.7, tcl=-.2, las=1)
plot(ns,abs(by.n['LR',]),type='l',ylim=range(abs(by.n)),xlab='n',ylab='distance to reference distribution')
lines(ns,abs(by.n['REML',]),lty=2)
legend("bottomleft",lty = 1:2,legend = c('LR test','IV test'),cex = 0.7,bty = "n",y.intersp=0.85,x.intersp=0.8)
par(op)
dev.off()



## plot routine for figures in the simulation section

save.dir <- '' # directory  with output of simulations
file.name <- '' # filename for tex output
ylab <- "rate to $\\hat{\\tau}^{2}_{MLE}$" # uncomment depending on sim
## ylab <- "rate to $T_{LR}$"
## ylab <- 'rate to reference distribution'

source('utils.R')
filelist <- dir(path=save.dir)
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]
filelist <- paste0(save.dir,filelist)
by.sim <- lapply(filelist, function(file) {
    load(file)
    list(by.n=by.n,ns=ns,m.exp=m.exp)
})
long <- lapply(by.sim, function(sim) {
    by.n <- sim$by.n
    ns <- sim$ns
    rates <- apply(by.n,1,function(gap)coef(lm(log(gap)~log(ns)))[2])
    settings <- c(m.exp=sim$m.exp)
    data.frame(estimator=names(rates),rate=unname(rates),t(settings))
})
long <- Reduce(rbind,long)
long <- aggregate(cbind(rate) ~ ., FUN=median, data=long)
print(length(filelist))
tikzDevice::tikz(file.name, width = 3, height = 2.5, pointsize = 8)
op <- par(mai=c(0.38, 0.48, 0.04, 0.03), mgp=c(1.35, 0.45, 0), cex.axis=0.60, cex.lab = 0.7, tcl=-.2, las=1)
plot(0, type = "n", xlim = range(long$m.exp), ylim = c(-2, .7), xlab = "$\\log_n m$", ylab = ylab)
estimators <- c("REML", "$S^2$", "DL", "HS", "HE", "PM", "SJ")
abline(h=-1)
lty <- c(1, 2, 3, 3, 4, 4, 4)
for (i in seq_along(estimators)) {
  lines(rate ~ m.exp, data = subset(long, estimator == estimators[i]), lty = lty[i])
}
legend("bottomleft",lty = sort(unique(lty)),legend = sapply(split(estimators, lty), paste, collapse = ", "),cex = 0.7,bty = "n",y.intersp=0.85,x.intersp=0.8)
par(op)
dev.off()
