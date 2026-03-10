library("hadron")
library("christianesfunctions")
library("optparse")
library("splines")


if (TRUE) {
  # set option list
  option_list <- list(
    make_option(c("-m", "--mode"), type = "character", default = "DG",
                help = "mode: DG, DM, DM2 [default %default]"),
    make_option(c("-k", "--kernel"), type = "character", default = "aic",
                help = "kernel [default %default]"),
    make_option(c("-f", "--folder"), type = "character", default = "-1",
                help = "folder with results [default %default]"),
    make_option(c("-p", "--plotfolder"), type = "character", default = "-1",
                help = "folder where plots and tables are stored [default %default]"),
    make_option(c("-b", "--bmass"), type = "integer", default = "-1",
                help = "index for bmass [default %default]"),
    make_option(c("-e", "--error"), type = "character", default = "stat",
                help = "index for bmass [default %default]"),
    make_option(c("-a", "--massfolder"), type = "character", default = "-1",
                help = "location where data for mass is stored [default %default]"),
    make_option(c("-l", "--ensemble"), type = "character", default = "-1",
                help = "ensemble that is analysed [default %default]")

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

if(opt$plotfolder=="-1") opt$plotfolder <- opt$folder
if(opt$massfolder=="-1") opt$massfolder <- opt$folder


errstring <- opt$error
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

## read in data
if(file.exists(sprintf("%s/%s_SCI_sigma_m%d_%s.csv", opt$folder, opt$mode, opt$bmass, errstring))) {
    res <- read.table(sprintf("%s/%s_SCI_sigma_m%d_%s.csv", opt$folder, opt$mode, opt$bmass, errstring), header=T)
} else if (file.exists(sprintf("%s/%s_SCI_sigma_m%d_%s_smear.csv", opt$folder, opt$mode, opt$bmass, errstring))) {
    res <- read.table(sprintf("%s/%s_SCI_sigma_m%d_%s_smear.csv", opt$folder, opt$mode, opt$bmass, errstring), header=T)
} else stop("results table cannot be read")

res <- res[res$kernel==opt$kernel, ]
masses <- read.table(sprintf("%s/masseslimitsingle.csv", opt$massfolder), header=T)
masses <- masses[masses$ensemble==opt$ensemble, ]
resint <- data.frame(iz=c(), int=c(), dint=c())
bsint <- array(NA, dim=c(zmax+2, 1000))
boot <- readRDS(sprintf("%s/%s_SCI_sigma_%s_m%d_%s.RDS", opt$folder, opt$mode, opt$kernel, opt$bmass, errstring))
# meanboot <- apply(boot, c(1, 2), mean)
# sdboot <- apply(boot, c(1, 2), mean)
# meanmean <- array(res$DG[res$mass==opt$bmass], dim=c(3, 10))
# print((meanmean - meanboot)/sdboot)
# res[res$mass==opt$bmass, ]

pdf(sprintf("%s/%s_SIC_integral_m%d_%s.pdf", opt$plotfolder, opt$mode, opt$bmass, errstring), title="")

omegamax <- 0.5*(1-(masses$M_F[opt$bmass]/masses$M_I[opt$bmass])^2)
xval <- (c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9.5)/10 * omegamax)^2
upperbound <- (omegamax)^2
for(iz in 0:zmax) {
  ## set boundaries, integral points

  yval <- c(0, res$DG[res$mass==opt$bmass & res$iz==iz & res$kernel==opt$kernel])

  bsamples <- array(NA, dim=c(1000, 11))
  bsamples[, 1] <- 0
  bsamples[, 2:11] <- t(boot[iz+1, ,])
  ## take out negative values
  # bsamples <- (bsamples + abs(bsamples))/2
  dyval <- apply(bsamples, 2, sd)

  ## plot splines
  xseq <- seq(0, upperbound, length.out=1000)
  plotwitherror(x=xval, y=yval, dy=dyval, xlab="omega^2", ylab="DGamma/domega^2", main=paste("spline interpolation mass", opt$bmass, "Z", iz), xlim=range(xval)*c(1, 1.2))
  abline(h=0)
  spline <- interpSpline(xval, yval)
  prediction <- t(apply(bsamples, 1, FUN=function(y) {
    inter <- interpSpline(xval, y)
    return(predict(object=inter, x=xseq)$y)
  }))
  # apply(prediction, 1, function(y) lines(x=xseq, y=y, col="red"))
  sds <- apply(prediction, 2, sd)
  means <- apply(prediction, 2, mean)
  pcol1 <- col2rgb("blue", alpha=TRUE)/255
  pcol1[4] <- 0.2
  pcol1 <- rgb(red=pcol1[1],green=pcol1[2],blue=pcol1[3],alpha=pcol1[4])

  polygon(x=c(xseq, rev(xseq)), y=c(means - sds, rev(means + sds)), col = pcol1)
  lines(x=xseq, y=predict(object=spline, x=xseq)$y, col="red", lty=1)
  lines(x=xseq, y=apply(prediction, 2, min), col="blue", lty=2)
  lines(x=xseq, y=apply(prediction, 2, max), col="blue", lty=2)
  plotwitherror(x=xval, y=yval, dy=dyval, rep=T)
  legend("topright", legend=c("meas", "mean", "68%", "min/max"), col=c("black", "red", pcol1, "blue"), pch=c(1, NA, 22, NA), pt.bg=c(NA, NA, pcol1, NA), lty=c(NA, 1, NA, 2))


  ## perform integral
  meanintspline <- splineintegral(yval=yval, xval=xval, continue = T,
                                  replacelower=F, higherlimit = upperbound,
                                  lowerlimit=upperbound, replaceindex=0)
  bsintspline <- apply(X=bsamples, MARGIN=1, FUN=splineintegral, xval=xval, continue = T,
                       replacelower=F, higherlimit = upperbound,
                       lowerlimit=upperbound, replaceindex=0)

  resint <- rbind(resint,data.frame(iz=iz, int=mean(bsintspline), dint=sd(bsintspline)) )
  bsint[iz+1, ] <- bsintspline
}

## sum of all Z-components

## set boundaries, integral points

yval <- c(0)
for(theta in unique(res$theta)) {
  yval <- append(yval, sum(res$DG[res$mass==opt$bmass & res$theta==theta & res$kernel==opt$kernel]))
}
bsamples <- array(NA, dim=c(1000, 11))
bsamples[, 1] <- 0
bsamples[, 2:11] <- t(apply(boot, c(2, 3), sum))
## take out negative values
# bsamples <- (bsamples + abs(bsamples))/2
dyval <- apply(bsamples, 2, sd)

## plot splines
xseq <- seq(0, upperbound, length.out=1000)
plotwitherror(x=xval, y=yval, dy=dyval, xlab="omega^2", ylab="DGamma/domega^2", main=paste("spline interpolation mass", opt$bmass, "Z", zmax+1), xlim=range(xval)*c(1, 1.2))
abline(h=0)
spline <- interpSpline(xval, yval)
prediction <- t(apply(bsamples, 1, FUN=function(y) {
  inter <- interpSpline(xval, y)
  return(predict(object=inter, x=xseq)$y)
}))
# apply(prediction, 1, function(y) lines(x=xseq, y=y, col="red"))
sds <- apply(prediction, 2, sd)
means <- apply(prediction, 2, mean)
pcol1 <- col2rgb("blue", alpha=TRUE)/255
pcol1[4] <- 0.2
pcol1 <- rgb(red=pcol1[1],green=pcol1[2],blue=pcol1[3],alpha=pcol1[4])

polygon(x=c(xseq, rev(xseq)), y=c(means - sds, rev(means + sds)), col = pcol1)
lines(x=xseq, y=predict(object=spline, x=xseq)$y, col="red", lty=1)
lines(x=xseq, y=apply(prediction, 2, min), col="blue", lty=2)
lines(x=xseq, y=apply(prediction, 2, max), col="blue", lty=2)
plotwitherror(x=xval, y=yval, dy=dyval, rep=T)
legend("topright", legend=c("meas", "mean", "68%", "min/max"), col=c("black", "red", pcol1, "blue"), pch=c(1, NA, 22, NA), pt.bg=c(NA, NA, pcol1, NA), lty=c(NA, 1, NA, 2))


## perform integral
meanintspline <- splineintegral(yval=yval, xval=xval, continue = T,
                                replacelower=F, higherlimit = upperbound,
                                lowerlimit=upperbound, replaceindex=0)
bsintspline <- apply(X=bsamples, MARGIN=1, FUN=splineintegral, xval=xval, continue = T,
                     replacelower=F, higherlimit = upperbound,
                     lowerlimit=upperbound, replaceindex=0)

resint <- rbind(resint,data.frame(iz=zmax+1, int=mean(bsintspline), dint=sd(bsintspline)) )
bsint[zmax+2, ] <- bsintspline

## save results
write.table(resint, sprintf("%s/%s_SIC_integral_%s_m%d_%s.csv", opt$plotfolder, opt$mode, opt$kernel, opt$bmass, errstring), row.names=F)
saveRDS(bsint, sprintf("%s/%s_SIC_integral_%s_m%d_%s.RDS", opt$plotfolder, opt$mode, opt$kernel, opt$bmass, errstring))
