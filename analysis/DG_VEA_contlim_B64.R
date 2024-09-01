
source("/hiskp4/gross/heavymesons/helpscripts/functions_stability_plots.R")
source("/hiskp4/gross/heavymesons/helpscripts/calc_DGDq2.R")
library("hadron")
zlist <- c(3, 0, 1, 2)
errlist <- c("stat", "sys", "vol", "tot")
doplot <- T
fnlin <- function(par, x, boot.R, ...) par[1] + par[2] * x
fncon <- function(par, x, boot.R, ...) par[1]

fitfn <- fnlin
par.guess <- c(1, 1)
comment <- "linear"

savefolder <- "tables_fnfour_10"

pdf("plots/DG_VEA_contlim.pdf", title="")


for(kernel in c("sigmoid", "erf")) {
  for(channel in c("cd", "cs")) {
    
    B64 <- readRDS(sprintf("%s/DG_VEA_epslim_converted_B64_%s_%s.RDS", savefolder, channel, kernel))
    C80 <- readRDS(sprintf("%s/DG_VEA_epslim_converted_C80_%s_%s.RDS", savefolder, channel, kernel))
    D96 <- readRDS(sprintf("%s/DG_VEA_epslim_converted_D96_%s_%s.RDS", savefolder, channel, kernel))
    E112 <- readRDS(sprintf("%s/DG_VEA_epslim_converted_E112_%s_%s.RDS", savefolder, channel, kernel))
    
    B64table <- read.table(sprintf("%s/DG_VEA_epslim_converted_B64_%s_%s.csv", savefolder, channel, kernel), header=TRUE)
    C80table <- read.table(sprintf("%s/DG_VEA_epslim_converted_C80_%s_%s.csv", savefolder, channel, kernel), header=TRUE)
    D96table <- read.table(sprintf("%s/DG_VEA_epslim_converted_D96_%s_%s.csv", savefolder, channel, kernel), header=TRUE)
    E112table <- read.table(sprintf("%s/DG_VEA_epslim_converted_E112_%s_%s.csv", savefolder, channel, kernel), header=TRUE)
    
    B64$L <- rep(64, length(B64$Nt))
    C80$L <- rep(80, length(C80$Nt))
    D96$L <- rep(96, length(D96$Nt))
    E112$L <- rep(112, length(E112$Nt))
    
    
    enslist <- list(B64, C80, D96, E112)
    tablelist <- list(B64table, C80table, D96table, E112table)
    
    
    res <- data.frame(th=NA, errtype=NA, iz=NA, lim=NA, dlim=NA, q=NA)
    reslist <- list(th=c(), errtype=c(), iz=c(), fit=list(), q=c())
    
    
    pcol <- col2rgb("gray", alpha=TRUE)/255 
    pcol[4] <- 0.65
    pcol <- rgb(red=pcol[1],green=pcol[2],blue=pcol[3],alpha=pcol[4])
    
    confcounter <- 1
    for (th in c(1, 2, 3, 4, 5, 6, 7, 8, 9, 9.5)){
      for (errtype in errlist) {
        for (iz in zlist) {
          title <- sprintf("%s %s th%s q^2 %.3f GeV^2 iz %d errtype %s", channel, kernel, th, B64table$q[abs(B64table$th - th) < 1e-2][1]^2, iz, errtype)
          print(title)
          ## we include the error on a in the continuum limit
          bsamples <- array(NA, dim=c(1000, 4))
          y <- c()
          dy <- c()
          x <- c()
          dx <- c()
          for (index in seq(1, 4)) {
            if(!is.na(tablelist[[index]]$q[abs(tablelist[[index]]$th - th) < 1e-2][1]^2)) {
              mask <- abs(enslist[[index]]$th - th) < 1e-2 & enslist[[index]]$iz == iz & enslist[[index]]$errtype==errtype
              mask[is.na(mask)] <- F
              masktable <- abs(tablelist[[index]]$th - th) < 1e-2 & tablelist[[index]]$iz == iz & tablelist[[index]]$errtype==errtype
              masktable[is.na(masktable)] <- F
              y[index] <- enslist[[index]]$DGDq2gev[mask]
              dy[index] <- enslist[[index]]$dDGDq2gev[mask]
              try(bsamples[, index] <- enslist[[index]]$datgev[, mask])
              x[index] <- tablelist[[index]]$afm[masktable]
              dx[index] <- tablelist[[index]]$dafm[masktable]
            }
          }
          
          
          
          fitresult <- try(bootstrap.nlsfit(x=x^2, y=y, bsamples=bsamples, fn=fitfn, par.guess=par.guess))
          
          if(!inherits(fitresult, "try-error")) {
            if(length(par.guess)==1) {
              if (doplot) plotwitherror(x=x^2, y=y, dy=dy, dx=dx, xlab="a[fm]^2", ylab="24 pi^3 DGamma / Dq^2 [GeV^-3]", main=title, xlim=c(-0.001, 0.007))
              if (doplot) polygon(x=c(-0.002, 0.01, 0.01, -0.002), 
                                  y=c(rep(fitresult$t0[1] - fitresult$se[1], 2), rep(fitresult$t0[1] + fitresult$se[1], 2)), 
                                  col=pcol, border=NA)
              if (doplot) lines(x=c(-0.002, 0.01), y=rep(fitresult$t0[1], 2))
              if (doplot) plotwitherror(x=x^2, y=y, dy=dy, rep=TRUE)
            } else {
              if (doplot) plot(fitresult, xlab="a[fm]^2", ylab="24 pi^3 DGamma / Dq^2 [GeV^-3]", 
                               main=title, xlim=c(-0.001, 0.007), plot.range=2*c(-0.001, 0.007))
            }
            if (doplot) plotwitherror(x=0, y=fitresult$t0[1], dy=fitresult$se[1], rep=TRUE, col="red")
            res <- rbind(res, data.frame(th=th, iz=iz, errtype=errtype, 
                                         lim=fitresult$t0[1], dlim=fitresult$se[1], q=B64table$q[abs(B64table$th - th) < 1e-2][1]))
            reslist$th[confcounter] <- th
            reslist$iz[confcounter] <- iz
            reslist$errtype[confcounter] <- errtype
            reslist$fit[[confcounter]] <- fitresult
            reslist$q[confcounter] <- B64table$q[abs(B64table$th - th) < 1e-2][1]
          } else {
            if (doplot) try(plotwitherror(x=x, y=y, dy=dy, xlab="a[fm]^2", ylab="24 pi^3 DGamma / Dq^2 [GeV^-3]", main=title))
            res <- rbind(res, data.frame(th=th, iz=iz, errtype=errtype, 
                                         lim=NA, dlim=NA, q=B64table$q[abs(B64table$th - th) < 1e-2][1]))
            reslist$th[confcounter] <- th
            reslist$iz[confcounter] <- iz
            reslist$errtype[confcounter] <- errtype
            reslist$fit[[confcounter]] <- fitresult
            reslist$q[confcounter] <- B64table$q[abs(B64table$th - th) < 1e-2][1] 
          }
          confcounter <- confcounter + 1
        }
      }
    }
    
    res <- res[-1, ]
    write.table(res, sprintf("%s/DG_VEA_contlim_%s_%s_%s.csv", savefolder, channel, kernel, comment))
    saveRDS(reslist, sprintf("%s/DG_VEA_contlim_%s_%s_%s.RDS", savefolder, channel, kernel, comment))
    
    
    
    res <- read.table(sprintf("%s/DG_VEA_contlim_%s_%s_%s.csv", savefolder, channel, kernel, comment), header=TRUE)
    for (errtype in errlist) {
      plot(NA,
           main=paste("continuum limit errtype", errtype, channel, kernel),
           xlim = c(0, max(res$q^2)), ylim = c(0, max(res$lim[res$errtype==errtype] + res$dlim[res$errtype==errtype])),
           xlab=expression(paste(q^2, " [", GeV^2, "]", sep="")),
           ylab="")
      mtext(side=2, text=expression(paste("24", pi^3, frac(D*Gamma , D*q^2),  "[", GeV^-3,"]")), line=1.5)
      for (iz in zlist) {
        mask <- res$iz==iz & res$errtype==errtype
        plotwitherror(x=res$q[mask]^2, y=res$lim[mask], dy=res$dlim[mask], 
                      rep=TRUE, col=iz+1)
      }
      
      legend(x = "topleft", 
             legend=seq(0, 3), col=seq(1, 4), pch=rep(1, 4), title="Z")
      ## sum up Z0+Z1+Z2
      dqsum <- rep(0, 10)
      ddqsum2 <- rep(0, 10)
      for (iz in c(0, 1, 2)) {
        mask <- res$iz==iz & res$errtype==errtype
        dqsum <- dqsum + res$lim[mask] 
        ddqsum2 <- ddqsum2 + res$dlim[mask]^2
      }
      mask <- res$iz==0 & res$errtype==errtype
      plotwitherror(x=res$q[mask]^2, y=dqsum, dy=sqrt(ddqsum2),
                    main=paste("continuum limit errtype", errtype, "cd sigmoid"),
                    xlim = c(0, max(res$q^2)), ylim = c(0, max(res$lim[res$errtype==errtype] + res$dlim[res$errtype==errtype])),
                    xlab=expression(paste(q^2, " [", GeV^2, "]", sep="")),
                    ylab="", col=4)
      mtext(side=2, text=expression(paste("24", pi^3, frac(D*Gamma , D*q^2),  "[", GeV^-3,"]")), line=1.5)
      for (iz in c(0, 1, 2)) {
        mask <- res$iz==iz & res$errtype==errtype
        plotwitherror(x=res$q[mask]^2, y=res$lim[mask], dy=res$dlim[mask], 
                      rep=TRUE, col=iz+1)
      }
      legend(x = "topleft", 
             legend=c(0, 1, 2, "sum"), col=seq(1, 4), pch=rep(1, 4), title="Z")
    }
  }
}

dev.off()
