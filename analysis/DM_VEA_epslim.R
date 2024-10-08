source("/hiskp4/gross/heavymesons/helpscripts/functions_stability_plots.R")
source("/hiskp4/gross/heavymesons/helpscripts/calc_DGDq2.R")
library("hadron")

fneight <- function(par, x, boot.R, ...) par[1] + par[2] * x**2 + par[3] * x**4 + par[4] * x**6 + par[5] * x**8
fnsix <- function(par, x, boot.R, ...) par[1] + par[2] * x**2 + par[3] * x**4 + par[4] * x**6
fnfour <- function(par, x, boot.R, ...) par[1] + par[2] * x**2 + par[3] * x**4

maxAmin <- 1

tsnk <- c(56, 56, 65, 78, 91, 48, 48, 48)
Nt <- c(32, 32, 37, 46, 53, 24, 24, 24)

afm <- c(0.07957, 0.07957, 0.06821, 0.05692, 0.04891, 0.07957, 0.07957, 0.07957)
amds <-  c(0.8, 0.8, 0.684, 0.57, 0.49, 0.8, 0.8, 0.8)
ensembles <- c("cB211.07.64", "cB211.07.96", "cC211.06.80_600", "cD211.054.96", "cE211.044.112_300", "cB211.07.48_300", "cB211.07.48_400", "cB211.07.64_48_36")
nameshort <- c("B64", "B96", "C80", "D96", "E112", "B48_300", "B48_400", "B64_48_36")
savefolder <- "tables_fnfour_15"
savefolder <- "tables_fnfour_12"

par.guess <- rep(1, 3)
fitfn <- fnfour

for(kernel in c("sigmoid", "erf")) {
  for(channel in c("cd", "cs")) {
    
    files <- c(sprintf("th1/%s_new/outputDMDq2/DMDq2.bin", kernel),
               sprintf("th2/%s_new/outputDMDq2/DMDq2.bin", kernel),
               sprintf("th3/%s_new/outputDMDq2/DMDq2.bin", kernel),
               sprintf("th4/%s_new/outputDMDq2/DMDq2.bin", kernel),
               sprintf("th5/%s_new/outputDMDq2/DMDq2.bin", kernel),
               sprintf("th6/%s_new/outputDMDq2/DMDq2.bin", kernel),
               sprintf("th7/%s_new/outputDMDq2/DMDq2.bin", kernel),
               sprintf("th8/%s_new/outputDMDq2/DMDq2.bin", kernel),
               sprintf("th9/%s_new/outputDMDq2/DMDq2.bin", kernel),
               sprintf("th9.5/%s_new/outputDMDq2/DMDq2.bin", kernel))
    
    th <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 9.5)
    nerr <- rep(2, 10)
    
    for(ens_index in seq_along(ensembles)) {
      pdf(sprintf("plots/DM_VEA_epslim_%s_%s_%s.pdf", nameshort[ens_index], channel, kernel))
      ens <- ensembles[ens_index]
      
      tsnk <- rep(tsnk[ens_index], 10)
      Nt <- rep(Nt[ens_index], 10)
      
      savename <- sprintf("%s/DM_VEA_epslim_%s_%s_%s", savefolder, nameshort[ens_index], channel, kernel)
      
      par.guess <- rep(1, 5)
      
      
      determineDGDq2_all(resultpath = sprintf("/hiskp4/gross/heavymesons/data/%s/%s/", channel, ens), filenames = files, 
                         tsnk = tsnk, Nt = Nt, th = th, nerr = nerr, amin = maxAmin, savename = savename,
                         fitfn = fitfn, par.guess = par.guess, errors=c("stat", "sys", "vol", "tot"),
                         volumetable=sprintf("tables/volume_interpolations/DM_tryB64factor_%s_%s.csv", channel, kernel),
                         neps=12, zlist=c(0, 1, 2, 3, 4), maxz=5, NDG=5)
      
      a <- afm[ens_index]
      agev <- afm[ens_index] / 0.1973269804 # fm / hbarc = GeV^-1
      da <- 0.00013
      am_H <- amds[ens_index] # dimensionless
      conversionfactor <- (am_H/agev)^5
      
      savename <- sprintf("%s/DM_VEA_epslim_%s_%s_%s", savefolder, nameshort[ens_index], channel, kernel)
      result <- read.table(paste0(savename, ".csv"), header=TRUE)
      
      result[, c("DGDq2gev", "dDGDq2gev")] <- result[, c("DGDq2", "dDGDq2")] * conversionfactor
      result$q <- result$w * am_H / agev
      result$afm <- rep(a, length(result$q))
      result$dafm <- rep(da, length(result$q))
      
      for (iz in c(4, 0, 1, 2, 3)) {
        for (err in c("stat", "sys", "vol", "tot")){
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
      savename <- sprintf("%s/DM_VEA_epslim_converted_%s_%s_%s", savefolder, nameshort[ens_index], channel, kernel)
      
      dat$DGDq2gev <- dat$DGDq2 * conversionfactor
      dat$dDGDq2gev <- dat$dDGDq2 * conversionfactor
      dat$datgev <- dat$dat * conversionfactor
      
      saveRDS(dat, paste0(savename, ".RDS"))
      write.table(result, paste0(savename, ".csv"))
      dev.off()
    }
  }
}
