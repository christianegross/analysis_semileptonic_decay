
source("/hiskp4/gross/heavymesons/helpscripts/functions_stability_plots.R")
source("/hiskp4/gross/heavymesons/helpscripts/calc_DGDq2.R")
source("/hiskp4/gross/heavymesons/helpscripts/functions_pull_factor.R")
library("hadron")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args)==1)
mode <- args[1]
stopifnot(mode=="DG" || mode=="DM" || mode=="DM2")

if(mode=="DG") {
  zlist <- c(0, 1, 2, 3)
  resbinname <- "DGammaDq2"
  NDG <- 4
  ylab <- "48 pi^3 / m_DS^5 DG / Dw^2 [GeV^-3]"
}

if(mode=="DM") {
  zlist <- c(0, 1, 2, 3, 4)
  resbinname <- "DMDq2"
  NDG <- 5
  ylab <- "96 pi^4 / m_DS^5 DM / Dw^2 [GeV^-2]"
}

if(mode=="DM2") {
  zlist <- c(0, 1, 2, 3, 4, 5)
  resbinname <- "DM2Dq2"
  NDG <- 6
  ylab <- "960 pi^5 / m_DS^5 DM2 / Dw^2 [GeV^-1]"
}
errlist <- c("stat")
sumup <-FALSE
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
fnlin <- function(par, x, boot.R, ...) par[1] + par[2] * x
fncon <- function(par, x, boot.R, ...) par[1] + 0*x

fitfn <- fnlin
par.guess <- c(1, 1)
comment <- "linear"

res <- data.frame(th=c(), epsilon=c(), iz=c(), lim=c(), dlim=c(), sloipe=c(), dslope=c(), chi=c(), 
                  ratio64=c(), dratio64=c(), ratio5644=c(), dratio5644=c(), ratio4836=c(), dratio4836=c(),
                  volcorr=c(), dvolcorr=c())
listvolcorr <- list()

pcol <- col2rgb("gray", alpha=TRUE)/255 
pcol[4] <- 0.65
pcol <- rgb(red=pcol[1],green=pcol[2],blue=pcol[3],alpha=pcol[4])

confcounter <- 1
epsilons <- seq(1, 20)-1 ## go to C-style indexing for reading in


for(channel in c("cd", "cs")) {
  for(kernel in c("sigmoid", "erf")) {
    for(th in c(seq(1, 9), 9.5)) {
      listvolcorr <- list(fitresult=list(), interpolate=list())
      enslist <- list()
      enslist$B48       <- read_in_DGDq2(resultpath=sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.48_400/th%s/%s_new/output%s/", channel, as.character(th), kernel, resbinname), filename=sprintf("%s.bin", resbinname), NDG=NDG)
      enslist$B64_56_44 <- read_in_DGDq2(resultpath=sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.64/th%s/%s_new/output%s/", channel, as.character(th), kernel, resbinname), filename=sprintf("%s.bin", resbinname), NDG=NDG)
      enslist$B64_48_36 <- read_in_DGDq2(resultpath=sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.64_48_36/th%s/%s_new/output%s/", channel, as.character(th), kernel, resbinname), filename=sprintf("%s.bin", resbinname), NDG=NDG)
      enslist$B96       <- read_in_DGDq2(resultpath=sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.96/th%s/%s_new/output%s/", channel, as.character(th), kernel, resbinname), filename=sprintf("%s.bin", resbinname), NDG=NDG)
      enslist$L <- c(48, 64, 64, 96)
      pdf(sprintf("plots/volume_%s_%s_%s.pdf", mode, channel, kernel), title="")
      for (icomb in c(0)) {
        for(iz in zlist) {
          for (epsilon in epsilons){
            
            title <- sprintf("th%s iz %d epsilon %d", th, iz, epsilon)
            print(title)
            bsamples <- array(NA, dim=c(1000, length(enslist)-1))
            y <- c()
            dy <- c()
            x <- c()
            for (index in seq(1, length(enslist)-1)) {
              tmp <- extractdata(data=enslist[[index]], iz=iz, epsilons = epsilon, icomb=icomb)
              y[index] <- tmp$DGammamean
              dy[index] <- tmp$dDGamma
              bsamples[, index] <- tmp$boot
            }
            x <- enslist$L
            fitresult <- try(bootstrap.nlsfit(x=x, y=y, bsamples=bsamples, fn=fitfn, par.guess=par.guess))
            
            if(!inherits(x = fitresult, what="try-error")){
              plot(fitresult,
                   main=paste(channel, kernel, "diff. decay rate epsilon", epsilon, ", z =", iz, "th", th),
                   xlab="L", ylab=ylab, xaxt="n")
              plotwitherror(x=x[3:4], y=y[3:4], dy=dy[3:4], rep=T, pch=1)
              ratio64 <- bsamples[, 2] / bsamples[, 3]
              ratio5644 <- bsamples[, 3] / bsamples[, 4]
              ratio4836 <- bsamples[, 1] / bsamples[, 2]
              VCDE <- predict(fitresult, x = 68.5)
              listvolcorr[["interpolate"]][[paste0("icomb", icomb, "iz", iz, "epsilon", epsilon)]] <- VCDE
              listvolcorr[["fitresult"]][[paste0("icomb", icomb, "iz", iz, "epsilon", epsilon)]] <- fitresult
              newline <- data.frame(th=th, epsilon=epsilon, iz=iz,
                                    lim=fitresult$t0[1], dlim=fitresult$se[1],
                                    slope=fitresult$t0[2], dslope=fitresult$se[2],
                                    chi=fitresult$chisqr / fitresult$dof,
                                    ratio64 = mean(ratio64), dratio64=sd(ratio64),
                                    ratio5644 = mean(ratio5644), dratio5644=sd(ratio5644),
                                    ratio4836 = mean(ratio4836), dratio4836=sd(ratio4836),
                                    volcorr=VCDE$val, dvolcorr=VCDE$err)
              res <- rbind(res, newline)
            } else {
              plotwitherror(x=x[3:4], y=y[3:4], dy=dy[3:4],
                            main=paste(channel, kernel, "diff. decay rate epsilon", epsilon, ", z =", iz, "th", th),
                            xlab="L", ylab=ylab,
                            ylim=c(min(y-dy),
                                   max(y+dy)),
                            xlim=c(1/96, 1/48), xaxt="n")
            }
            axis(side=1, at=c(48, 64, 96), label=c("48", "64", "96"))
            legend(x="topleft", legend=c("56, 44", "48, 36", "pred."), col=c(1, 2, 4), pch=c(1, 1, 2), title=c("tsnk, tj"))
            plotwitherror(x=x[1:2], y=y[1:2], dy=dy[1:2],
                          rep=TRUE, col=2, pch=1, lwd=1.2)
            plotwitherror(x=VCDE$x, y=VCDE$val, dy=VCDE$err,
                          rep=TRUE, col=4, pch=1, lwd=1.2)
            
          }
        }
      }
      listvolcorr$neps <- length(epsilons)
      listvolcorr$epsilons <- enslist[[1]][[2]]$epsilons[epsilons+1]
      listvolcorr$Z <- zlist
      print(names(listvolcorr))
      saveRDS(object=listvolcorr, file=sprintf("/hiskp4/gross/heavymesons/analysis/tables/volume_interpolations/%s_vol_B_th%s_%s_%s.RDS", mode, as.character(th), channel, kernel))
      rm(enslist)
    }
  }
}


pdf(sprintf("plots/%s_volerror.pdf", mode), title="")
for(channel in c("cd", "cs")) {
  for(kernel in c("sigmoid", "erf")) {
    filenames <- rep(sprintf("%s.bin", resbinname), 10)
    resultpathB96 <- c(sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.96/th1/%s_new/output%s/", channel, kernel, resbinname),
                       sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.96/th2/%s_new/output%s/", channel, kernel, resbinname),
                       sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.96/th3/%s_new/output%s/", channel, kernel, resbinname),
                       sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.96/th4/%s_new/output%s/", channel, kernel, resbinname),
                       sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.96/th5/%s_new/output%s/", channel, kernel, resbinname),
                       sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.96/th6/%s_new/output%s/", channel, kernel, resbinname),
                       sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.96/th7/%s_new/output%s/", channel, kernel, resbinname),
                       sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.96/th8/%s_new/output%s/", channel, kernel, resbinname),
                       sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.96/th9/%s_new/output%s/", channel, kernel, resbinname),
                       sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.96/th9.5/%s_new/output%s/", channel, kernel, resbinname))
    
    resultpathB64 <- c(sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.64/th1/%s_new/output%s/", channel, kernel, resbinname),
                       sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.64/th2/%s_new/output%s/", channel, kernel, resbinname),
                       sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.64/th3/%s_new/output%s/", channel, kernel, resbinname),
                       sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.64/th4/%s_new/output%s/", channel, kernel, resbinname),
                       sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.64/th5/%s_new/output%s/", channel, kernel, resbinname),
                       sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.64/th6/%s_new/output%s/", channel, kernel, resbinname),
                       sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.64/th7/%s_new/output%s/", channel, kernel, resbinname),
                       sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.64/th8/%s_new/output%s/", channel, kernel, resbinname),
                       sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.64/th9/%s_new/output%s/", channel, kernel, resbinname),
                       sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.64/th9.5/%s_new/output%s/", channel, kernel, resbinname))
    
    filenamesinterpolated <- c(sprintf("%s_vol_B_th1_%s_%s.RDS", mode, channel, kernel),
                               sprintf("%s_vol_B_th2_%s_%s.RDS", mode, channel, kernel),
                               sprintf("%s_vol_B_th3_%s_%s.RDS", mode, channel, kernel),
                               sprintf("%s_vol_B_th4_%s_%s.RDS", mode, channel, kernel),
                               sprintf("%s_vol_B_th5_%s_%s.RDS", mode, channel, kernel),
                               sprintf("%s_vol_B_th6_%s_%s.RDS", mode, channel, kernel),
                               sprintf("%s_vol_B_th7_%s_%s.RDS", mode, channel, kernel),
                               sprintf("%s_vol_B_th8_%s_%s.RDS", mode, channel, kernel),
                               sprintf("%s_vol_B_th9_%s_%s.RDS", mode, channel, kernel),
                               sprintf("%s_vol_B_th9.5_%s_%s.RDS", mode, channel, kernel))
    
    resultpathinterpolated <- rep("/hiskp4/gross/heavymesons/analysis/tables/volume_interpolations/", 10)
    
    nerr <- rep(2, 10)
    L96 <- rep(96, 10)
    L64 <- rep(64, 10)
    Linter <- rep(68.5, 10)
    th <- c(seq(1, 9), 9.5)
    
    
    dat1 <- determinesyserrfinitevolumeextrapolated(filenames1 = filenames, 
                                                    resultpath1 = resultpathB96, 
                                                    nerr1 = nerr, L1 = L96, 
                                                    resultpath2 = resultpathinterpolated, 
                                                    filenames2 = filenamesinterpolated, nerr2 = nerr, L2 = Linter, th = th, 
                                                    savename = sprintf("/hiskp4/gross/heavymesons/analysis/tables/volume_interpolations/%s_trypullfactor_%s_%s", mode, channel, kernel), 
                                                    mode=mode)
    print(dat1)
    
    
    dat2 <- determinesyserrfinitevolume(filenames1 = filenames, 
                                        resultpath1 = resultpathB96, 
                                        nerr1 = nerr, L1 = L96, 
                                        resultpath2 = resultpathB64, 
                                        filenames2 = filenames, nerr2 = nerr, L2 = L64, th = th, 
                                        savename = sprintf("/hiskp4/gross/heavymesons/analysis/tables/volume_interpolations/%s_tryB64factor_%s_%s", mode, channel, kernel), 
                                        mode=mode)
    print(dat2)
    
    for(iz in zlist) {
      plot(x=seq(1, length(dat1$P[dat1$iz==iz])), y=dat1$P[dat1$iz==iz],
           pch=1, col=1, xlab="index", ylab="P", main=paste(channel, kernel, "Z", iz),
           ylim=range(na.omit(dat2$P[dat2$iz==iz], dat1$P[dat1$iz==iz])))
      points(x=seq(1, length(dat2$P[dat2$iz==iz])), y=dat2$P[dat2$iz==iz], col=2, pch=2)
      legend(x="topleft", legend=c("inter", "B64"), horiz=T, col=c(1, 2), pch=c(1, 2))
    }
  }
}
