#~ ---
#~ title: "Ratio method decay rate"
#~ format: pdf
#~ editor: source
#~ author: "Christiane"
#~ date: "2026-01-15"
#~ toc: true
#~ number-sections: true
#~ execute: 
#~   echo: false
#~   warning: false
#~   message: false
#~ ---

#~ ```{r setup, include=FALSE}
library("hadron")
library("christianesfunctions")
#~ getresultsatone<- function(x, na.rm=T) {
#~   fitfn <- x$fn
#~   data.frame(x=1, y=fitfn(x$t0[1:length(x$par.guess)], 1), 
#~              dy=sd(apply(x$t[, 1:length(x$par.guess)], MARGIN=1, fitfn, x=1), na.rm=na.rm))
#~ }
#~ mbs <- 5.36691

erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
#~ ```

#~ We perform the limits in the order HLT $\to$ continuum  $\to$ smearing $\to$ integrate $\to$ m. 
#~ Exchanging the order of continuum and sigma can in principle be done, it is important to consider the Finite Volume Effects before doing the smearing limit. However we are neglecting Finite Volume Effects so far.

#~ Here we consider if the continuum limit can be stabilised by using the ratio method.

#~ As an example we consider theta 5, Z0, sigma 0 and the sigmoid kernel.

#~ We consider first doing the continuum limit in $Gamma$ and then computing $$\frac{\Gamma_i m_{i+1}^5}{\Gamma_{i+1}m_i^5}$$.
#~ We also consider first computing the ratios and then doing the continuum limit of the ratios.

#~ ```{r}
args <- commandArgs(trailingOnly = TRUE)
th <- as.character(args[1])
iz <- as.numeric(args[2])
isigma <- as.numeric(args[3])

pdf(sprintf("Documents/heavymesons/bdecay/continuum/plotsratio/ratio_%s_z%d_sigma%d.pdf", th, iz, isigma), title="")
ensembles <- c("cB211.07.64", "cC211.06.80", "cD211.054.96")
afm <- c(0.07957, 0.06821, 0.05692)
agev <- afm / 0.1973269804
datboot <- array(NA, dim=c(6, 3, 1001))
massboot <- array(NA, dim=c(1001, 6, 3))
masscontlimit <- array(NA, dim=c(1001, 6))
tablemasses <- read.table("~/Documents/heavymesons/bdecay/continuum/tables/masseslimit.csv", header=T)
masscontlimit[1, ] <- tablemasses$M_I
tablemasses <- readRDS("~/Documents/heavymesons/bdecay/continuum/tables/masseslimit.RDS")
masscontlimit[2:1001, ] <- tablemasses[, 1:6]
for(mass in 1:6) {
  for(ens in ensembles){
    ## Gamma 
    tmp <- read.table(sprintf("~/Documents/heavymesons/data/bdecay/btoc_spec_s/%s/%s/plots_backcomparison_3_sqrtsum/chosen_BJnu_m%d_%s_ik0.csv", ens, th, mass, th), header=T)
    datboot[mass, which(ens==ensembles), 1] <- tmp$rho[tmp$iz==iz & tmp$ieps==isigma]
    datboot[mass, which(ens==ensembles), 2:1001] <- readRDS(sprintf("~/Documents/heavymesons/data/bdecay/btoc_spec_s/%s/%s/plots_backcomparison_3_sqrtsum/chosen_BJnu_m%s_%s_ik0_boots_err_stat_HLT.RDS", ens, th, mass, th))[iz+1, isigma+1, ]
    
    ## mass
    tmp <- read_in_Y(filename="Y.bin", resultpath=sprintf("~/Documents/heavymesons/data/bdecay/btoc_spec_s/%s/%s/m%d%s/outputY/", ens, th, mass, ifelse(ens=="cC211.06.80", "_sigmaphys", "")))[[2]]
    tmpboot <- tmp$mH$bsamples / agev[which(ens==ensembles)]
    massboot[, mass, which(ens==ensembles)] <- tmpboot[c(1002, 1:1000)]
  }
}
# dat

#~ First continuum limit, then ratio

fitfn <- function(par, x, boot.r, ...) par[1] + par[2]*x
contlimit <- array(NA, dim=c(1001, 6))
for(mass in 1:6) {
  fit <- bootstrap.nlsfit(x=afm^2, y=datboot[mass, , 1], bsamples = t(datboot[mass, ,2:1001]), fn = fitfn, par.guess=c(1, 1))
  plot(fit, xlab="a^2 [fm^2]", ylab="Gamma/Gammabar", main=paste("mass", mass, "first cont limit"), plot.range=c(0, max(afm^2)), 
       xlim=c(0, max(afm^2)), ylim=range(fit$y-fit$dy, fit$y+fit$dy, fit$t0[1]+c(1, -1)*fit$se[1]))
  plotwitherror(x=0, y=fit$t0[1], dy=fit$se[1], col="red", pch=0, rep=T)
  contlimit[1, mass] <- fit$t0[1]
  contlimit[2:1001, mass] <- fit$t[,1]
  pulllin <- (fit$t0[1]-fit$y[which.min(fit$x)])/fit$se[1]
  dsyslin <- abs(fit$t0[1]-fit$y[which.min(fit$x)])*erf(abs(pulllin)/sqrt(2))
  contlimit[2:1001, mass] <- contlimit[2:1001, mass] + rnorm(1000, mean=0, sd=dsyslin)
}
ratioaftercontlim <- contlimit[, 1:5] / contlimit[, 2:6] * (masscontlimit[, 2:6]/ masscontlimit[, 1:5])^5
ratioaftercontlim[1, ]
# plotwitherror(x=1:5, y=ratioaftercontlim[1, ], dy=apply(ratioaftercontlim[2:1001, ], 2, sd), xlab="index", ylab="ratio(Gamma/m^5)")


#~ First ratio, then continuum limit

fitfn <- function(par, x, boot.r, ...) par[1] + par[2]*x
ratio <- array(NA, dim=c(1001, 5, 3))
contlimafterratio <- array(NA, dim=c(1001, 5))
for(ens in ensembles){
  ratio[, , which(ens==ensembles)] <- t(datboot[1:5, which(ens==ensembles), ]) / t(datboot[2:6, which(ens==ensembles), ]) * (massboot[, 2:6, which(ens==ensembles)]/ massboot[, 1:5, which(ens==ensembles)])^5
}
for(index in 1:5) {
  fit <- bootstrap.nlsfit(x=afm^2, y=ratio[1, index, ], bsamples = ratio[2:1001, index, ], fn = fitfn, par.guess=c(1, 1))
  plot(fit, xlab="a^2 [fm^2]", ylab="ratio Gamma/m^5", main=paste("index", index, "first ratio"), plot.range=c(0, max(afm^2)), 
       xlim=c(0, max(afm^2)), ylim=range(fit$y-fit$dy, fit$y+fit$dy, fit$t0[1]+c(1, -1)*fit$se[1]))
  plotwitherror(x=0, y=fit$t0[1], dy=fit$se[1], col="red", pch=0, rep=T)
  contlimafterratio[1, index] <- fit$t0[1]
  contlimafterratio[2:1001, index] <- fit$t[,1]
  pulllin <- (fit$t0[1]-fit$y[which.min(fit$x)])/fit$se[1]
  dsyslin <- abs(fit$t0[1]-fit$y[which.min(fit$x)])*erf(abs(pulllin)/sqrt(2))
  contlimafterratio[2:1001, index] <- contlimafterratio[2:1001, index] + rnorm(1000, mean=0, sd=dsyslin)
}
# plotwitherror(x=1:5, y=contlimafterratio[1, ], dy=apply(contlimafterratio[2:1001, ], 2, sd), xlab="index", ylab="ratio(Gamma/m^5)")

#~ Compare


plotwitherror(x=1:5, y=contlimafterratio[1, ], dy=apply(contlimafterratio[2:1001, ], 2, sd), xlab="index", ylab="$$\\frac{\\Gamma_i m_{i+1}^5}{\\Gamma_{i+1}m_i^5}$$", xlim=c(1, 6))
plotwitherror(x=1:5+0.1, y=ratioaftercontlim[1, ], dy=apply(ratioaftercontlim[2:1001, ], 2, sd), col="red", pch=0, rep=T)
legend("right", legend=c("cont first", "ratio first"), pch=c(0, 1), col=c("red", "black"))
