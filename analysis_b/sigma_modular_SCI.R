library("hadron")

library("optparse")
library("christianesfunctions")


if (TRUE) {
  # set option list
  option_list <- list(
    make_option(c("-m", "--mode"), type = "character", default = "DG",
                help = "mode: DG, DM, DM2 [default %default]"),
    make_option(c("-f", "--folder"), type = "character", default = "-1",
                help = "folder with results [default %default]"),
    make_option(c("-p", "--plotfolder"), type = "character", default = "-1",
                help = "folder where plots and tables are stored [default %default]"),
    make_option(c("-l", "--ensemble"), type = "character", default = "-1",
                help = "ensemble to look at [default %default]"),
    make_option(c("-b", "--bmass"), type = "integer", default = "-1",
                help = "index for bmass [default %default]"),
    make_option(c("-e", "--error"), type = "character", default = "stat",
                help = "index for bmass [default %default]"),
    make_option(c("-s", "--subfolder"), type = "character", default = "plots_backcomparison_3_sqrtsum",
                help = "subfolder where results are stored [default %default]")
    
  )
  parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
  args <- parse_args(parser, positional_arguments = 0)
  opt <- args$options
}


erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

if(opt$mode=="DG") {
  modefolder <- "DGammaDq2"
  zmax <- 2
} else if(opt$mode=="DM") {
  modefolder <- "DMDq2"
  zmax <- 3
} else if(opt$mode=="DM2") {
  modefolder <- "DM2Dq2"
  zmax <- 4
}


afm <- c(0.07957, 0.06821, 0.05692)
ensembles <- c("cB211.07.64", "cC211.06.80", "cD211.054.96")
stopifnot(opt$ensemble %in% ensembles)
afm <- afm[which(opt$ensemble == ensembles)]

if(opt$plotfolder=="-1") opt$plotfolder <- opt$folder


fnfourcomb <- function (par, x, boot.r, maskfn, ...) {
  return (par[1]  +  par[2] * x^2 * maskfn + par[3] * x^4 * maskfn + par[4] * x^2 * (!maskfn) + par[5] * x^4 * (!maskfn))
} 
fnsixcomb <- function (par, x, boot.r, maskfn, ...) {
  return (par[1]  +  par[2] * x^4 * maskfn + par[3] * x^6 * maskfn + par[4] * x^4 * (!maskfn) + par[5] * x^6 * (!maskfn))
}       

fitfnz01 <- function(pars, x, boot.r, ...) pars[1] + x^2*pars[2] + x^4*pars[3]
fitfnz23 <- function(pars, x, boot.r, ...) pars[1] + x^2*pars[2] + x^4*pars[3]

# maskfn <- c(rep(T, length(dattmpsigmoid$eps)), rep(F, length(dattmperf$eps)))
#         myfit <- try(bootstrap.nlsfit(fn=fnfourcomb, x=c(dattmpsigmoid$eps, dattmperf$eps), 
#                                       y=c(dattmpsigmoid$DGDq2, dattmperf$DGDq2), bs=bscomb, 
#                                       par.guess=c(1, 1, 1,1, 1), mask=maskfit, 
#                                       maskfn=maskfn[maskfit], na.rm=T))

errstring <- opt$error
thetas <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 9.5)
# thetas <- c(1)

erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

res <- data.frame(ensemble=c(), mass=c(), iz=c(), theta=c(), DG=c(), dDG=c(), pull=c(), dsys=c(), kernel=c(), weight=c())


bootssigmoid <- array(NA, dim=c(3, 10, 1000))
bootssigmoid_sys <- array(NA, dim=c(3, 10, 1000))
bootserf <- array(NA, dim=c(3, 10, 1000))
bootserf_sys <- array(NA, dim=c(3, 10, 1000))
bootscomb <- array(NA, dim=c(3, 10, 1000))
bootscomb_sys <- array(NA, dim=c(3, 10, 1000))
bootsaic <- array(NA, dim=c(3, 10, 1000))
bootsaic_sys <- array(NA, dim=c(3, 10, 1000))

pdf(sprintf("%s/%s_SCI_sigma_m%d_%s_smear.pdf", opt$plotfolder, opt$mode, opt$bmass, errstring), title="")
for(theta in thetas) {
  if(file.exists(sprintf("%s/th%s/%s/%schosen_BJnu_m%d_th%s_ik0.csv", opt$folder, as.character(theta), opt$subfolder, ifelse(opt$mode=="DG", "", paste0("_", opt$mode)), opt$bmass, as.character(theta)))) {
    tablesigmoid <- read.table(sprintf("%s/th%s/%s/%schosen_BJnu_m%d_th%s_ik0.csv", opt$folder, as.character(theta), opt$subfolder, ifelse(opt$mode=="DG", "", paste0("_", opt$mode)), opt$bmass, as.character(theta)), header=T)
    tableerf <- read.table(sprintf("%s/th%s/%s/%schosen_erf_BJnu_m%d_th%s_ik0.csv", opt$folder, as.character(theta), opt$subfolder, ifelse(opt$mode=="DG", "", paste0("_", opt$mode)), opt$bmass, as.character(theta)), header=T)
  } else stop("result table cannot be read")
  
  tablesigmoid$sigma <- tablesigmoid$eps / afm * 0.1973269804
  tableerf$sigma <- tableerf$eps / afm * 0.1973269804
  
  bootsigmoid <- readRDS(sprintf("%s/th%s/%s/%schosen_BJnu_m%d_th%s_ik0_boots_err_%s.RDS", opt$folder, as.character(theta), opt$subfolder, ifelse(opt$mode=="DG", "", paste0("_", opt$mode)), opt$bmass, as.character(theta), errstring))
  booterf <- readRDS(sprintf("%s/th%s/%s/%schosen_erf_BJnu_m%d_th%s_ik0_boots_err_%s.RDS", opt$folder, as.character(theta), opt$subfolder, ifelse(opt$mode=="DG", "", paste0("_", opt$mode)), opt$bmass, as.character(theta), errstring))
  
  for(iz in 0:2) {
    fitfncomb <- ifelse(iz >= 2, fnsixcomb, fnfourcomb)
    fitfn <- ifelse(iz >= 2, fitfnz23, fitfnz01)
    
    fitsigmoid <- bootstrap.nlsfit(fn=fitfn, par.guess=c(1, 1, 1), y=tablesigmoid$rho[tablesigmoid$iz==iz], x=tablesigmoid$sigma[tablesigmoid$iz==iz], bsamples=t(bootsigmoid[iz+1,, ]), success.infos = 1:4, na.rm=T)
    fiterf <- bootstrap.nlsfit(fn=fitfn, par.guess=c(1, 1, 1), y=tableerf$rho[tableerf$iz==iz], x=tableerf$sigma[tableerf$iz==iz]/2.5, bsamples=t(booterf[iz+1,, ]), success.infos = 1:4, na.rm=T)
    
    maskfn <- c(rep(T, 10), rep(F, 8))
    bscomb <- array(c(t(bootsigmoid[iz+1,, ]), t(booterf[iz+1,, ])), dim=c(1000, 18))
    
    fitcomb <- bootstrap.nlsfit(fn=fitfncomb, x=c(tablesigmoid$sigma[tablesigmoid$iz==iz], tableerf$sigma[tableerf$iz==iz]/2.5), 
                                y=c(tablesigmoid$rho[tablesigmoid$iz==iz], tableerf$rho[tableerf$iz==iz]), bs=bscomb, 
                                par.guess=c(1, 1, 1,1, 1), 
                                maskfn=maskfn, na.rm=T)
    
    
    listfits <- list(fitsigmoid, fiterf, fitcomb)
    weights <- sapply(listfits, function(x){
      exp(-0.5*(x$chisqr + 2*length(x$par.guess) - length(x$x)))
    })
    weights <- weights/sum(weights)
    
    average <- sum(sapply(listfits, function(x) x$t0[1])*weights)
    averageboot <- apply(sapply(listfits, function(x) x$t[, 1]), 1, function(x) sum(x*weights))
    averagese <- sd(averageboot)
    
    plotcombined(fitresult=fitcomb, xlim=c(-0.003, max(fitcomb$x)), 
                 ylim=range(fitcomb$y-fitcomb$dy, fitcomb$y+fitcomb$dy)*c(0.9, 1.2), 
                 cols=c("black", "black"),
                 xlab="sigma [GeV]", ylab="DG/Dq2", pchs=c(1, 1), polygon=F,
                 main=paste(opt$ensemble, "mass", opt$bmass, "theta", as.character(theta), "iz", iz))
    plot(fitsigmoid, col.band="red", col="red", col.line="red", rep=T,
         plot.range=c(0, max(fitsigmoid$x)), opacity.band=0.0)
    plot(fiterf, col.band="blue", col="blue", col.line="blue", rep=T,
         plot.range=c(0, max(fiterf$x)), opacity.band=0.0)
    plotwitherror(x=-0.001, y=fitsigmoid$t0[1], dy=fitsigmoid$se[1], col="red", pch=20, rep=T, cex=2)
    plotwitherror(x=-0.002, y=fiterf$t0[1], dy=fitsigmoid$se[1], col="blue", pch=20, rep=T, cex=2)
    plotwitherror(x=-0.003, y=fitcomb$t0[1], dy=fitcomb$se[1], col="black", pch=20, rep=T, cex=2)
    plotwitherror(x=0, y=average, dy=averagese, col="darkgreen", pch=2, lwd=3, rep=T)
    legend(x="top", legend=paste(c("sigmoid", "erf", "comb", "res"), c(format(weights, digits=2), "")), 
           col=c("red", "blue", "black", "darkgreen"), pch=c(1,1,1,2))
    
    ## set pull factors and determine systematic errors
    ## if only one fit, use pull factor to estimatesys error
    ## if multiple values, determine sys error by calculating total error first and assuming errors add in quadrature
    pullsigmoid <- (fitsigmoid$t0[1]-fitsigmoid$y[1])/fitsigmoid$se[1]
    dsyssigmoid <- abs(fitsigmoid$t0[1]-fitsigmoid$y[1])*erf(abs(pullsigmoid)/sqrt(2))
    res <- rbind(res, data.frame(mass=opt$bmass, iz=iz, theta=theta, DG=fitsigmoid$t0[1], dDG=fitsigmoid$se[1], pull=pullsigmoid, dsys=dsyssigmoid, kernel="sigmoid", weight=weights[1]))
    
    pullerf <- (fiterf$t0[1]-fiterf$y[1])/fiterf$se[1]
    dsyserf <- abs(fiterf$t0[1]-fiterf$y[1])*erf(abs(pullerf)/sqrt(2))
    res <- rbind(res, data.frame(mass=opt$bmass, iz=iz, theta=theta, DG=fiterf$t0[1], dDG=fiterf$se[1], pull=pullerf, dsys=dsyserf, kernel="erf", weight=weights[2]))
    
    pullcomb <- c((average-fiterf$y[which.min(fiterf$x)])/averagese, (average-fitsigmoid$y[which.min(fitsigmoid$x)])/fitcomb$se)
    maxindex <- which.max(abs(pullcomb))
    pullcomb <- pullcomb[maxindex]
    dsyscomb <- abs(fitcomb$t0[1]-c(fiterf$y[which.min(fiterf$x)], fitsigmoid$y[which.min(fitsigmoid$x)])[maxindex])*erf(abs(pullcomb)/sqrt(2))
    res <- rbind(res, data.frame(mass=opt$bmass, iz=iz, theta=theta, DG=fitcomb$t0[1], dDG=fitcomb$se[1], pull=pullcomb, dsys=dsyscomb, kernel="comb", weight=weights[3]))
    
    ## for AIC, pull factor can come from sigmoid and erf
    pullaic <- c((average-fiterf$y[which.min(fiterf$x)])/averagese, (average-fitsigmoid$y[which.min(fitsigmoid$x)])/averagese)
    pullaic <- pullaic[which.max(abs(pullaic))]
    dtotaic <- sqrt(sum(sapply(listfits, function(x) x$se[1]^2)*weights) + sum(sapply(listfits, function(x) (x$t0[1]-average)^2)*weights))
    dsysaic <- sqrt(dtotaic^2-averagese^2)
    res <- rbind(res, data.frame(mass=opt$bmass, iz=iz, theta=theta, DG=average, dDG=averagese, pull=pullaic, dsys=dsysaic, kernel="aic", weight=1))
    
    bootssigmoid[iz+1, which(theta==thetas), ] <- fitsigmoid$t[, 1]
    bootssigmoid_sys[iz+1, which(theta==thetas), ] <- fitsigmoid$t[, 1] + rnorm(1000, 0, dsyssigmoid)
    bootserf[iz+1, which(theta==thetas), ] <- fiterf$t[, 1]
    bootserf_sys[iz+1, which(theta==thetas), ] <- fiterf$t[, 1] + rnorm(1000, 0, dsyserf)
    bootscomb[iz+1, which(theta==thetas), ] <- fitcomb$t[, 1]
    bootscomb_sys[iz+1, which(theta==thetas), ] <- fitcomb$t[, 1] + rnorm(1000, 0, dsyscomb)
    bootsaic[iz+1, which(theta==thetas), ] <- averageboot
    bootsaic_sys[iz+1, which(theta==thetas), ] <- averageboot + rnorm(1000, 0, dsyserf)
    
  }
  
}
saveRDS(bootssigmoid, file=sprintf("%s/%s_SCI_sigma_%s_m%d_%s.RDS", opt$plotfolder, opt$mode, "sigmoid", opt$bmass, errstring))
saveRDS(bootssigmoid_sys, file=sprintf("%s/%s_SCI_sigma_%s_m%d_%s_smear.RDS", opt$plotfolder, opt$mode, "sigmoid", opt$bmass, errstring))
saveRDS(bootserf, file=sprintf("%s/%s_SCI_sigma_%s_m%d_%s.RDS", opt$plotfolder, opt$mode, "erf", opt$bmass, errstring))
saveRDS(bootserf_sys, file=sprintf("%s/%s_SCI_sigma_%s_m%d_%s_smear.RDS", opt$plotfolder, opt$mode, "erf", opt$bmass, errstring))
saveRDS(bootscomb, file=sprintf("%s/%s_SCI_sigma_%s_m%d_%s.RDS", opt$plotfolder, opt$mode, "comb", opt$bmass, errstring))
saveRDS(bootscomb_sys, file=sprintf("%s/%s_SCI_sigma_%s_m%d_%s_smear.RDS", opt$plotfolder, opt$mode, "comb", opt$bmass, errstring))
saveRDS(bootsaic, file=sprintf("%s/%s_SCI_sigma_%s_m%d_%s.RDS", opt$plotfolder, opt$mode, "aic", opt$bmass, errstring))
saveRDS(bootsaic_sys, file=sprintf("%s/%s_SCI_sigma_%s_m%d_%s_smear.RDS", opt$plotfolder, opt$mode, "aic", opt$bmass, errstring))
write.table(res, sprintf("%s/%s_SCI_sigma_m%d_%s_smear.csv", opt$plotfolder, opt$mode, opt$bmass, errstring), row.names=F)
