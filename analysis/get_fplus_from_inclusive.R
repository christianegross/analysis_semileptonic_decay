library("hadron")
library("christianesfunctions")
library("optparse")


if (TRUE) {
  # set option list
  option_list <- list(
    make_option(c("-c", "--channel"), type = "character", default = "cs",
                help = "decay channel [default %default]"),
    make_option(c("-e", "--ensemble"), type = "character", default = "cB211.07.64",
                help = "ensemble [default %default]"),
    make_option(c("-t", "--momentum"), type = "character", default = "th2",
                help = "momentum [default %default]"),
    make_option(c("-f", "--folder"), type = "character", default = "-1",
                help = "folder with results [default %default]"),
    make_option(c("-p", "--plotfolder"), type = "character", default = "-1",
                help = "folder where plots and tables are stored [default %default]"),
    make_option(c("-b", "--bmass"), type = "integer", default = "-1",
                help = "index for bmass [default %default]"),
    make_option(c("-i", "--icomb"), type = "integer", default = "0",
                help = "icomb: (1, FULL), (2, VPAR), (3, APAR), (4, VPERP), (5, APERP) [default %default]"),
    make_option(c("--t1"), type = "integer", default = "10",
                help = "lower limit plateau [default %default]"),
    make_option(c("--t2"), type = "integer", default = "20",
                help = "upper limit plateau [default %default]"),
    make_option(c("-r", "--rfg"), type = "double", default = "0.5",
                help = "M_final/M_initial [default %default]")
    
  )
  parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
  args <- parse_args(parser, positional_arguments = 0)
  opt <- args$options
}
if(opt$folder=="-1") {
  opt$folder <- sprintf("~/Documents/heavymesons/data/%s/%s/%s", opt$channel, opt$ensemble, opt$momentum)
}
if(opt$plotfolder=="-1") {
  opt$plotfolder <- opt$folder
}


bmassaddon <- ifelse(opt$bmass==-1, "", sprintf("_bmass_%d", opt$bmass))


pdf(sprintf("%s/4Pt_inclusive%s_%s.pdf", opt$plotfolder, bmassaddon, opt$momentum), title="", height=11.7, width=16.6)


setwd("~/Documents/heavymesons/scripts/analysis/")

## icomb: (1, FULL), (2, VPAR), (3, APAR), (4, VPERP), (5, APERP)
icomb <- opt$icomb
th <- opt$momentum
rfg <- opt$rfg
resultpath <- opt$folder
t1 <- opt$t1+1
t2 <- opt$t2+1
stopifnot(t1 < t2)
stopifnot(0 <= t1)


## read into list
dat <- read_in_Y(resultpath=resultpath, filename="Y.bin")[[2]]
mH <- dat$mH$bsamples[dat$mH$n+2]
wmin <- sqrt(rfg^2+dat$metadata$w^2)

stopifnot(t2 <= dat$metadata$tj)


## make more handy arrays: Ymean[iy, icomb, t], ..., Yboot[iy, icomb, t, iboot]
Yar <- makeYarray(dat)


## construct correlators according to eq. 104 of 2504.06063v2

# Yplusmean <- (Yar$Ymean[2,icomb , ] + (1-wmin)^2/(dat$metadata$w^2)*Yar$Ymean[3,icomb , ] - 2*(1-wmin)/dat$metadata$w*Yar$Ymean[4,icomb , ])*(2*pi/mH )
Yplusmwob <- (Yar$Ymwob[2,icomb , ] + (1-wmin)^2/(dat$metadata$w^2)*Yar$Ymwob[3,icomb , ] - 2*(1-wmin)/dat$metadata$w*Yar$Ymwob[4,icomb , ])*(2*pi/mH )
## bootstrap samples: make correct array and then take the individual mass of each sample
Yplusdata <- t(Yar$Ydata[2,icomb,, ] + (1-wmin)^2/(dat$metadata$w^2)*Yar$Ydata[3,icomb,, ] - 2*(1-wmin)/dat$metadata$w*Yar$Ydata[4,icomb,, ])*(2*pi)
dim(Yplusdata)
length(dat$mH$bsamples[1:dat$mH$n])
Yplusdata <- sweep(Yplusdata, 1, dat$mH$bsamples[1:dat$mH$n], FUN = '/')
Yplusmean <- apply(Yplusdata, 2, mean)
Yplussd <- apply(Yplusdata, 2, sd)
Yplusbias <- Yplusmean - Yplusmwob

## effective residue as in fig 1 of 2504.06063v2
## Why is there a factor of M_H in the exponent, but not in the factor in addition to those in the paper?

# effresmean <- wmin/2/pi * exp(wmin*((0:dat$metadata$tj))*mH ) * Yplusmean
effresmwob <- wmin/2/pi * exp(wmin*((0:dat$metadata$tj))*mH ) * Yplusmwob
effresdata <- wmin/2/pi * exp(wmin*outer(dat$mH$bsamples[1:dat$mH$n], (0:dat$metadata$tj)))  * Yplusdata
effresmean <- apply(effresdata, 2, mean)
effressd <- apply(effresdata, 2, sd)
effresbias <- effresmean - effresmwob


# plot

plotwitherror(x=0:dat$metadata$tj, y=Yplusmwob, dy=Yplussd, log="y", main="Y^+", xlab="t/a")
plotwitherror(x=0:dat$metadata$tj, y=effresmwob, dy=effressd, main="[f^+]^2", xlab="t/a", ylab="f^+^2")
plotwitherror(x=0:dat$metadata$tj, y=effresmwob, dy=effressd, main=sprintf("[f^+]^2 mass %d %s", opt$bmass, opt$momentum), 
              xlab="t/a", ylab="f^+^2", ylim=mean(effresmwob[6:dat$metadata$tj][effresmwob>1e-2])+c(-1, 1)*sd(effresmwob[6:dat$metadata$tj][effresmwob>1e-2]))
grid()
# plotwitherror(x=0:dat$metadata$tj, y=abs(Yar$Ymwob[1,icomb,]), dy=Yar$Ysd[1,icomb,], log="y", main="Y1", xlab="t/a")
# plotwitherror(x=0:dat$metadata$tj, y=abs(Yar$Ymwob[2,icomb,]), dy=Yar$Ysd[2,icomb,], log="y", main="Y2", xlab="t/a")
# plotwitherror(x=0:dat$metadata$tj, y=abs(Yar$Ymwob[3,icomb,]), dy=Yar$Ysd[3,icomb,], log="y", main="Y3", xlab="t/a")
# plotwitherror(x=0:dat$metadata$tj, y=abs(Yar$Ymwob[4,icomb,]), dy=Yar$Ysd[4,icomb,], log="y", main="Y4", xlab="t/a")
# plotwitherror(x=0:dat$metadata$tj, y=abs(Yar$Ymwob[5,icomb,]), dy=Yar$Ysd[5,icomb,], log="y", main="Y5", xlab="t/a")

## fit effective residue to a constant in region t1 to t2
fnconst <- function(par, x, boot.r, ...) par[1] + 0*x

fit.result <- bootstrap.nlsfit(fn = fnconst, par.guess=c(1), y=effresmwob, x=0:dat$metadata$tj, bsamples = effresdata, mask=t1:t2)

plot(fit.result, 
     main=sprintf("[f^+]^2=%s mass %d %s", tex.catwitherror(fit.result$t0[1], fit.result$se[1], with.dollar=F, with.cdot=F, digits=2), opt$bmass, opt$momentum), 
     xlab="t/a", ylab="f^+^2",
     ylim=range(effresmwob[t1:t2])*c(0.85, 1.15))
grid()
abline(v=1:20*5, col="gray", lty=3)

## eq. A.37 of DsXlnu notes by Nazario
## leave out factor Gammabar=mH^5*G_F^2*S_EW/(48*pi^4) to have better control at varying B masses

dgdw2mean <- 2*pi/wmin*dat$metadata$w^3*fit.result$t0[1]
dgdw2boot <- 2*pi/wmin*dat$metadata$w^3*fit.result$t[, 1]
dgdw2sd <- sd(dgdw2boot)

tex.catwitherror(dgdw2mean, dgdw2sd)

res <- list(fit=fit.result, Yplus=list(mean=Yplusmean, mwob=Yplusmwob, bias=Yplusbias, sd=Yplussd, data=Yplusdata),
            w=dat$metadata$w, wmin=wmin, mH=mH, rfg=rfg, icomb=icomb, t1=t1, t2=t2,
            dgdw2 = list(data=dgdw2boot, mean=mean(dgdw2boot), mwob=dgdw2mean, bias=mean(dgdw2boot)-dgdw2mean, sd=dgdw2sd))

saveRDS(res, file=sprintf("%s/4Pt_inclusive%s_%s.RDS", opt$plotfolder, bmassaddon, opt$momentum))


dev.off()



