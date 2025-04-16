source("/hiskp4/gross/heavymesons/helpscripts/functions_stability_plots.R")
source("/hiskp4/gross/heavymesons/helpscripts/calc_DGDq2.R")
library("hadron")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args)==1 || length(args)==2)
mode <- args[1]
stopifnot(mode=="DG" || mode=="DM" || mode=="DM2")
if(length(args)==2) stopifnot(args[2]=="su")

if(mode=="DG") {
    zlist <- c(0, 1, 2, 3)
    maxz <- 4
    NDG <- 4
    output <- "DGammaDq2"
}

if(mode=="DM") {
    zlist <- c(0, 1, 2, 3, 4)
    maxz <- 5
    NDG <- 5
    output <- "DMDq2"
}

if(mode=="DM2") {
    zlist <- c(0, 1, 2, 3, 4, 5)
    maxz <- 6
    NDG <- 6
    output <- "DM2Dq2"
}

docdcs <- T
dosu <- F
if(length(args)==2) {
    docdcs <- F
    dosu <- T
}



fneight <- function(par, x, boot.R, ...) par[1] + par[2] * x**2 + par[3] * x**4 + par[4] * x**6 + par[5] * x**8
fnsix <- function(par, x, boot.R, ...) par[1] + par[2] * x**2 + par[3] * x**4 + par[4] * x**6
fnfour <- function(par, x, boot.R, ...) par[1] + par[2] * x**2 + par[3] * x**4

maxAmin <- 1

tsnkall <- c(56, 56, 65, 78, 91, 48, 48, 48)
Ntall <- c(32, 32, 37, 46, 53, 24, 24, 24)

afm <- c(0.07957, 0.07957, 0.06821, 0.05692, 0.04891, 0.07957, 0.07957, 0.07957)
amds <-  c(0.8, 0.8, 0.684, 0.57, 0.49, 0.8, 0.8, 0.8)
ensembles <- c("cB211.07.64", "cB211.07.96", "cC211.06.80_600", "cD211.054.96", "cE211.044.112_300", "cB211.07.48_300", "cB211.07.48_400", "cB211.07.64_48_36")
nameshort <- c("B64", "B96", "C80", "D96", "E112", "B48_300", "B48_400", "B64_48_36")
savefolder <- "tables_fnfour_15"
savefolder <- "tables_fnfour_20"

errlist <- c("stat", "sys", "vol", "tot")

par.guess <- rep(1, 3)
fitfn <- fnfour

if(docdcs) {
for(kernel in c("sigmoid", "erf")) {
  for(channel in c("cd", "cs")) {
    
    files <- c(sprintf("th1/%s_new/output%s/%s.bin", kernel, output, output),
               sprintf("th2/%s_new/output%s/%s.bin", kernel, output, output),
               sprintf("th3/%s_new/output%s/%s.bin", kernel, output, output),
               sprintf("th4/%s_new/output%s/%s.bin", kernel, output, output),
               sprintf("th5/%s_new/output%s/%s.bin", kernel, output, output),
               sprintf("th6/%s_new/output%s/%s.bin", kernel, output, output),
               sprintf("th7/%s_new/output%s/%s.bin", kernel, output, output),
               sprintf("th8/%s_new/output%s/%s.bin", kernel, output, output),
               sprintf("th9/%s_new/output%s/%s.bin", kernel, output, output),
               sprintf("th9.5/%s_new/output%s/%s.bin", kernel, output, output))
    
    th <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 9.5)
    nerr <- rep(2, 10)
    
    for(ens_index in seq_along(ensembles)) {
      pdf(sprintf("plots/%s_VEA_epslim_%s_%s_%s.pdf", mode, nameshort[ens_index], channel, kernel))
      ens <- ensembles[ens_index]
      
      tsnk <- rep(tsnkall[ens_index], 10)
      Nt <- rep(Ntall[ens_index], 10)
      
      savename <- sprintf("%s/%s_VEA_epslim_%s_%s_%s", savefolder, mode, nameshort[ens_index], channel, kernel)
      
      
      determineDGDq2_all(resultpath = sprintf("/hiskp4/gross/heavymesons/data/%s/%s/", channel, ens), filenames = files, 
                         tsnk = tsnk, Nt = Nt, th = th, nerr = nerr, amin = maxAmin, savename = savename,
                         fitfn = fitfn, par.guess = par.guess, errors=errlist,
                         volumetable=sprintf("tables/volume_interpolations/%s_tryB64factor_%s_%s.csv", mode, channel, kernel),
                         neps=20, zlist=zlist, maxz=maxz, NDG=NDG)
      
      a <- afm[ens_index]
      agev <- afm[ens_index] / 0.1973269804 # fm / hbarc = GeV^-1
      da <- 0.00013
      am_H <- amds[ens_index] # dimensionless
      conversionfactor <- 1
      
      
      savename <- sprintf("%s/%s_VEA_epslim_%s_%s_%s", savefolder, mode, nameshort[ens_index], channel, kernel)
      result <- read.table(paste0(savename, ".csv"), header=TRUE)
      
      result[, c("DGDq2gev", "dDGDq2gev")] <- result[, c("DGDq2", "dDGDq2")] * conversionfactor
      result$q <- result$w * am_H / agev
      result$afm <- rep(a, length(result$q))
      result$dafm <- rep(da, length(result$q))
      
      for (iz in zlist) {
        for (err in errlist){
          mask <- result$iz==iz & result$icomb==0 & result$errtype==err
          try(plotwitherror(x=result$q[mask]^2, y=result$DGDq2gev[mask], dy=result$dDGDq2gev[mask],
                            main=paste("differential decay rate for tsink=56, tins=44, sigmoid, cd, ", err, "Z", iz), 
                            xlab="q^2 [GeV^2]", ylab="96 pi^4 DM / Dq^2 [GeV^-3]",
                            ylim=c(0, max(result$DGDq2gev[mask] + result$dDGDq2gev[mask])), xlim=c(0, max(result$q[mask]^2))))
          
          try(plotwitherror(x=result$q[mask]^2, y=result$dDGDq2gev
                            [mask],
                            main=paste("differential decay rate for tsink=NA, tins=44, sigmoid, cd, ", err, "Z", iz), 
                            ylab="96 pi^4 d(DM / Dq^2) [GeV^-3]", xlim=c(0, max(result$q[mask]^2))))
          
        }
      }
      
      dat <- readRDS(paste0(savename, ".RDS"))
      savename <- sprintf("%s/%s_VEA_epslim_converted_%s_%s_%s", savefolder, mode, nameshort[ens_index], channel, kernel)
      
      write.table(result, paste0(savename, ".csv"))
      dev.off()
    }
  }
}
}


## su does not have all ensembles and different momenta
tsnkall <- c(56, 78)
Ntall <- c(32, 46)

afm <- c(0.07957, 0.05692)
amds <-  c(0.8, 0.57)
ensembles <- c("cB211.07.64", "cD211.054.96")
nameshort <- c("B64", "D96")
savefolder <- "tables_fnfour_20"

errlist <- c("stat", "sys")

par.guess <- rep(1, 3)
fitfn <- fnfour

if(dosu) {
for(kernel in c("sigmoid", "erf")) {
  for(channel in c("su")) {
    
    filesall <- c(sprintf("th2_su/%s_new/output%s/%s.bin", kernel, output, output),
               sprintf("th4_su/%s_new/output%s/%s.bin", kernel, output, output),
               sprintf("th6_su/%s_new/output%s/%s.bin", kernel, output, output),
               sprintf("th8_su/%s_new/output%s/%s.bin", kernel, output, output),
               sprintf("th9.5_su/%s_new/output%s/%s.bin", kernel, output, output), 
               sprintf("th10_su/%s_new/output%s/%s.bin", kernel, output, output),
               sprintf("th1/%s_new/output%s/%s.bin", kernel, output, output))
    
    thall <- c("th2_su", "th4_su", "th6_su", "th8_su", "th9.5_su", "th10_su", "th1")
    nerrall <- rep(2, 7)
    
    
    for(ens_index in seq_along(ensembles)) {
      if(nameshort[ens_index]=="B64") {
          files <- filesall[c(1:5, 7)]
          th <- thall[c(1:5, 7)]
          nerr <- nerrall[c(1:5, 7)]
      } else if (nameshort[ens_index]=="D96") {
          files <- filesall[c(1:6)]
          th <- thall[c(1:6)]
          nerr <- nerrall[c(1:6)]
      }
      pdf(sprintf("plots/%s_VEA_epslim_%s_%s_%s.pdf", mode, nameshort[ens_index], channel, kernel))
      ens <- ensembles[ens_index]
      
      tsnk <- rep(tsnkall[ens_index], 6)
      Nt <- rep(Ntall[ens_index], 6)
      
      savename <- sprintf("%s/%s_VEA_epslim_%s_%s_%s", savefolder, mode, nameshort[ens_index], channel, kernel)
      
      
      determineDGDq2_all(resultpath = sprintf("/hiskp4/gross/heavymesons/data/%s/%s/", channel, ens), filenames = files, 
                         tsnk = tsnk, Nt = Nt, th = th, nerr = nerr, amin = maxAmin, savename = savename,
                         fitfn = fitfn, par.guess = par.guess, errors=errlist,
                         neps=20, zlist=zlist, maxz=maxz, NDG=NDG)
      
      a <- afm[ens_index]
      agev <- afm[ens_index] / 0.1973269804 # fm / hbarc = GeV^-1
      da <- 0.00013
      am_H <- amds[ens_index] # dimensionless
      conversionfactor <- 1
      
      
      savename <- sprintf("%s/%s_VEA_epslim_%s_%s_%s", savefolder, mode, nameshort[ens_index], channel, kernel)
      result <- read.table(paste0(savename, ".csv"), header=TRUE)
      
      result[, c("DGDq2gev", "dDGDq2gev")] <- result[, c("DGDq2", "dDGDq2")] * conversionfactor
      result$q <- result$w * am_H / agev
      result$afm <- rep(a, length(result$q))
      result$dafm <- rep(da, length(result$q))
      
      print(result)
      
      for (err in errlist){
        mask <- result$icomb==0 & result$errtype==err
        try(plotwitherror(x=result$w[mask]^2+0.00005*result$iz[mask], y=result$DGDq2[mask], dy=result$dDGDq2[mask], 
                          xlab="w^2", ylab=sprintf("%s/Dw^2", mode),
                          ylim=range(result$DGDq2[mask] + result$dDGDq2[mask], result$DGDq2[mask] - result$dDGDq2[mask]), xlim=c(0, max(result$w[mask]^2)), 
                          col=result$iz[mask]+1, pch=result$iz[mask]+1))
        try(legend(x="topleft", legend=zlist, col=zlist+1, pch=zlist+1))
        try(plotwitherror(x=result$w[mask]^2+0.00005*result$iz[mask], y=result$dDGDq2[mask], 
                          xlab="w^2", ylab=sprintf("%s/Dw^2", mode),
                          ylim=c(0, max(result$dDGDq2[mask])), xlim=c(0, max(result$w[mask]^2)), 
                          col=result$iz[mask]+1, pch=result$iz[mask]+1))
        try(legend(x="topleft", legend=zlist, col=zlist+1, pch=zlist+1))
        
        
      }
      
      dat <- readRDS(paste0(savename, ".RDS"))
      savename <- sprintf("%s/%s_VEA_epslim_converted_%s_%s_%s", savefolder, mode, nameshort[ens_index], channel, kernel)
      
      write.table(result, paste0(savename, ".csv"))
      dev.off()
    }
  }
}
}
