library("hadron")

library("optparse")


if (TRUE) {
  # set option list
  option_list <- list(
    make_option(c("-m", "--mode"), type = "character", default = "DG",
                help = "mode: DG, DM, DM2 [default %default]"),
    make_option(c("-k", "--kernel"), type = "character", default = "aic",
                help = "kernel [default %default]"),
    make_option(c("-i", "--integraltype"), type = "character", default = "averagewithsyserr",
                help = "type of integral [default %default]"),
    make_option(c("-f", "--folder"), type = "character", default = "-1",
                help = "folder with results [default %default]"),
    make_option(c("-s", "--subfolder"), type = "character", default = "analyse/tables",
                help = "subfolder with results after ensemble specification [default %default]"),
    make_option(c("-p", "--plotfolder"), type = "character", default = "-1",
                help = "folder where plots and tables are stored [default %default]"),
    make_option(c("-b", "--bmass"), type = "integer", default = "-1",
                help = "index for bmass [default %default]"),
    make_option(c("-e", "--error"), type = "character", default = "stat",
                help = "errors included in data [default %default]")

  )
  parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
  args <- parse_args(parser, positional_arguments = 0)
  opt <- args$options
}

errstring <- opt$error
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

integraltypes <- c("spline", "trapezoidal", "simpson", "average", "averagewithsyserr")
stopifnot(opt$integraltype %in% integraltypes)

## the results for integraltypes spline, trapezoid and simpson are stored together, so have to be extracted differently
additionalint <- opt$integraltype %in% c("spline", "trapezoidal", "simpson")
if(additionalint & grepl("int", errstring)) stop("there can be no systematic error for one integral type")


if(additionalint) {
    B64 <- readRDS(sprintf("%s/B64/%s/%s_SIC_integral_%s_m%d_%s_all.RDS", opt$folder, opt$subfolder, opt$mode, opt$kernel, opt$bmass, errstring))[which(opt$integraltype==integraltypes),,]
    nboot <- dim(B64)[2]
    
    C80 <- readRDS(sprintf("%s/C80/%s/%s_SIC_integral_%s_m%d_%s_all.RDS", opt$folder, opt$subfolder, opt$mode, opt$kernel, opt$bmass, errstring))[which(opt$integraltype==integraltypes),,]
    stopifnot(nboot == dim(C80)[2])
    
    D96 <- readRDS(sprintf("%s/D96/%s/%s_SIC_integral_%s_m%d_%s_all.RDS", opt$folder, opt$subfolder, opt$mode, opt$kernel, opt$bmass, errstring))[which(opt$integraltype==integraltypes),,]
    stopifnot(nboot == dim(D96)[2])
} else {

    B64 <- readRDS(sprintf("%s/B64/%s/%s_SIC_integral_%s_m%d_%s.RDS", opt$folder, opt$subfolder, opt$mode, opt$kernel, opt$bmass, errstring))
    nboot <- dim(B64)[2]
    
    C80 <- readRDS(sprintf("%s/C80/%s/%s_SIC_integral_%s_m%d_%s.RDS", opt$folder, opt$subfolder, opt$mode, opt$kernel, opt$bmass, errstring))
    stopifnot(nboot == dim(C80)[2])
    
    D96 <- readRDS(sprintf("%s/D96/%s/%s_SIC_integral_%s_m%d_%s.RDS", opt$folder, opt$subfolder, opt$mode, opt$kernel, opt$bmass, errstring))
    stopifnot(nboot == dim(D96)[2])
}

## table with mean values only has type average with stat and sys error, no separate line for sys.
tableintegral <- opt$integraltype
if(opt$integraltype=="averagewithsyserr") tableintegral <- "average"

if(file.exists(sprintf("%s/B64/%s/%s_SIC_integral_%s_m%d_%s.csv", opt$folder, opt$subfolder, opt$mode, opt$kernel, opt$bmass, errstring))) {
    B64_table <- read.table(sprintf("%s/B64/%s/%s_SIC_integral_%s_m%d_%s.csv", opt$folder, opt$subfolder, opt$mode, opt$kernel, opt$bmass, errstring), header=T)
    C80_table <- read.table(sprintf("%s/C80/%s/%s_SIC_integral_%s_m%d_%s.csv", opt$folder, opt$subfolder, opt$mode, opt$kernel, opt$bmass, errstring), header=T)
    D96_table <- read.table(sprintf("%s/D96/%s/%s_SIC_integral_%s_m%d_%s.csv", opt$folder, opt$subfolder, opt$mode, opt$kernel, opt$bmass, errstring), header=T)
} else if (file.exists(sprintf("%s/B64/%s/%s_SIC_integral_%s_m%d_%s_int.csv", opt$folder, opt$subfolder, opt$mode, opt$kernel, opt$bmass, errstring))) {
    B64_table <- read.table(sprintf("%s/B64/%s/%s_SIC_integral_%s_m%d_%s_int.csv", opt$folder, opt$subfolder, opt$mode, opt$kernel, opt$bmass, errstring), header=T)
    C80_table <- read.table(sprintf("%s/C80/%s/%s_SIC_integral_%s_m%d_%s_int.csv", opt$folder, opt$subfolder, opt$mode, opt$kernel, opt$bmass, errstring), header=T)
    D96_table <- read.table(sprintf("%s/D96/%s/%s_SIC_integral_%s_m%d_%s_int.csv", opt$folder, opt$subfolder, opt$mode, opt$kernel, opt$bmass, errstring), header=T)
} else stop("results table cannot be read")



B64_table <- B64_table[B64_table$type==tableintegral, ]
C80_table <- C80_table[C80_table$type==tableintegral, ]
D96_table <- D96_table[D96_table$type==tableintegral, ]

reslist <- list()
reslistlinear <- list()
reslistsys <- list()
reslistlinearsys <- list()
restable <- data.frame(iz=c(), DG=c(), dDG=c(), dsys=c(), dtot=c(), pull=c(), cutoff=c(), extrapolation=c())


fitfn1 <- function(pars, x, boot.r, ...) pars[1] + x*pars[2]
fitfn2 <- function(pars, x, boot.r, ...) pars[1] + x*0

pdf(sprintf("%s/%s_SIC_contlimit_%s_%s_m%d_%s.pdf", opt$plotfolder, opt$mode, opt$kernel, opt$integraltype, opt$bmass, errstring))

afm <- c(0.07957, 0.06821, 0.05692) ## superseeded, from 2206.15084
afm <- c(0.07948, 0.06819, 0.05685) ## from 2411.08852
dafm <- c(1.1, 1.4, 0.9)*1e-4
afmbootsamples <- parametric.bootstrap(1000, afm, dafm, 1234)^2
for(iz in 0:(zmax+1)) {

  boots <- array(c(B64[iz+1, ], C80[iz+1, ], D96[iz+1, ], afmbootsamples), dim=c(nboot, 6))
  y <- c(B64_table$int[B64_table$iz==iz], C80_table$int[C80_table$iz==iz], D96_table$int[D96_table$iz==iz])

  # perform fits and assign AIC weights
  fit1 <- bootstrap.nlsfit(x=afm^2, y=y, bsamples=boots, fn=fitfn1, par.guess=c(1, 1), success.infos=1:4)
  weights <- c(exp(-1*(fit1$chisqr + 2 * 2 - length(fit1$x))/2))
  fit2 <- try(bootstrap.nlsfit(x=afm^2, y=y, bsamples=boots, fn=fitfn2, par.guess=c(1), mask=c(F, T, T), success.infos=1:4))
  if(inherits(fit2, "try-error")){
    weights[2] <- 0
    fit2 <- list(t0=c(0), se=c(0), t=array(rep(0, nboot), dim=c(nboot, 1)), failed=T)
  } else {
    weights[2] <- c(exp(-1*(fit2$chisqr + 2 * 1 - length(fit1$x))/2))
    fit2$failed=F
  }

  ## AIC
  weights <- weights/sum(weights)
  # print(weights)
  average <- sum(c(fit1$t0[1], fit2$t0[1])*weights)
  averageboot <- fit1$t[, 1]*weights[1] + fit2$t[, 1]*weights[2]
  averagese <- sd(averageboot, na.rm=T)

  ## plot with Histogramm of Bootstrap samples at the left side
  layout.matrix <- matrix(c(1, 2, 2, 2), nrow = 1, ncol = 4)
  layout(mat = layout.matrix,
         heights = 1, # Heights of the two rows
         widths = 1) # Widths of the two columns

  mai <- par("mai")
  par(mai=c(mai[1], mai[2], mai[3], 0))

  myhist <- hist(averageboot, plot = FALSE, breaks=50)
  #~     myhist <- hist(c(fit1$t[, 1], fit2$t[, 1]), plot = FALSE, breaks=50)
  nbars <- length(myhist$counts)
  midbreaks <- (myhist$breaks[2:(nbars+1)] + myhist$breaks[1:nbars]) / 2
  plot(x=myhist$counts, y=midbreaks,
       xlim=c(max(myhist$counts), min(myhist$counts)),type="l",
       ylim=range(myhist$breaks, fit1$t0[1] + c(-1, 1)*fit1$se[1], fit1$y+fit1$dy, fit1$y-fit1$dy),
       xaxt="n", xlab="", ylab="DGDq2")

  par(mai=c(mai[1], 0, mai[3], mai[2]))

  plot(fit1,
       xlab="a^2[fm^2]", ylab="", main=paste("Z", iz),
       plot.range=c(0, 0.07957^2), xlim=c(0, 0.07957^2),
       ylim=range(myhist$breaks, fit1$t0[1] + c(-1, 1)*fit1$se[1], fit1$y+fit1$dy, fit1$y-fit1$dy),
       yaxt="n", col="black", col.band="black", col.line="black", opacity.band=weights[1])
  if(!fit2$failed) {
    plot(fit2, rep=T,
         plot.range=c(0, 0.06821^2), xlim=c(0, 0.06821^2),
         yaxt="n", col="blue", col.band="blue", col.line="blue", opacity.band=weights[2])
  }

  plotwitherror(x=0, y=average, dy=averagese, col="red", pch=21, rep=T)
  plotwitherror(x=0.0002, y=fit1$t0[1], dy=fit1$se[1], col="magenta", pch=21, rep=T)
  plotwitherror(x=0.0002, y=fit2$t0[1], dy=fit2$se[1], col="magenta", pch=21, rep=T)
  legend("top", legend=c(paste("linear, w=", format(weights[1], digits=2)), paste("constant, w=", format(weights[2], digits=2)), "result", "single results"),
         col=c("black", "blue", "red", "magenta"), pch=c(1, 1, 21, 21), lty=c(1, 1, NA, NA))
  axis(side=4)
  par(mai=mai)

  pull <- (average-fit1$y[3])/averagese
  ## change this to be like in paper
  dsys <- abs(average-fit1$y[3])*erf(abs(pull)/sqrt(2))
  dtot <- sqrt(sum(c(fit1$se[1], fit2$se[1])^2*weights) + sum(c(fit1$t0[1]-average, fit2$t0[1]-average)^2*weights))
  dsys <- sqrt(dtot^2-averagese^2)

  res <- list(mean=average, sd=averagese, boot=averageboot, pull=0, cutoff=(average-fit1$y[3])/average,
              dsys=0, dtot=averagese, iz=iz)
  reslist[[paste0("iz", iz)]] <- res
  res <- list(mean=average, sd=averagese, boot=averageboot + rnorm(nboot, 0, dsys), pull=pull, cutoff=(average-fit1$y[3])/average,
              dsys=dsys, dtot=dtot, iz=iz)
  reslistsys[[paste0("iz", iz)]] <- res

  pulllin <- (fit1$t0[1]-fit1$y[3])/fit1$se[1]
  dsyslin <- abs(fit1$t0[1]-fit1$y[3])*erf(abs(pulllin)/sqrt(2))

  reslin <- list(mean=fit1$t0[1], sd=fit1$se[1], boot=fit1$t[, 1], pull=0, cutoff=(fit1$t0[1]-fit1$y[3])/fit1$t0[1],
                 dsys=0, dtot=fit1$se[1], iz=iz)
  reslistlinear[[paste0("iz", iz)]] <- reslin
  reslin <- list(mean=fit1$t0[1], sd=fit1$se[1], boot=fit1$t[, 1] + rnorm(nboot, 0, dsyslin), pull=pulllin, cutoff=(fit1$t0[1]-fit1$y[3])/fit1$t0[1],
                 dsys=dsyslin, dtot=sqrt(fit1$se[1]^2+dsyslin^2), iz=iz, fit=fit1)
  reslistlinearsys[[paste0("iz", iz)]] <- reslin

  #~     res <- data.frame(iz=c(), isigma=c(), DG=c(), dDG=c(), dsys=c(), dtot=c(), pull=c(), cutoff=c(), extrapolation=c())

  restable <- rbind(restable, data.frame(iz=iz, DG=average, dDG=averagese, dsys=dsys, dtot=dtot, pull=pull, cutoff=(average-fit1$y[3])/average, extrapolation="aic"))
  restable <- rbind(restable, data.frame(iz=iz, DG=fit1$t0[1], dDG=fit1$se[1], dsys=dsyslin, dtot=sqrt(fit1$se[1]^2+dsyslin^2), pull=pulllin, cutoff=(fit1$t0[1]-fit1$y[3])/fit1$t0[1], extrapolation="linear"))
}

saveRDS(object=reslist,          file=sprintf("%s/%s_SIC_contlimit_%s_%s_m%d_%s_AIC.RDS",         opt$plotfolder, opt$mode, opt$kernel, opt$integraltype, opt$bmass, errstring))
saveRDS(object=reslistlinear,    file=sprintf("%s/%s_SIC_contlimit_%s_%s_m%d_%s_linear.RDS",      opt$plotfolder, opt$mode, opt$kernel, opt$integraltype, opt$bmass, errstring))
saveRDS(object=reslistsys,       file=sprintf("%s/%s_SIC_contlimit_%s_%s_m%d_%s_cont_AIC.RDS",    opt$plotfolder, opt$mode, opt$kernel, opt$integraltype, opt$bmass, errstring))
saveRDS(object=reslistlinearsys, file=sprintf("%s/%s_SIC_contlimit_%s_%s_m%d_%s_cont_linear.RDS", opt$plotfolder, opt$mode, opt$kernel, opt$integraltype, opt$bmass, errstring))
write.table(restable,            file=sprintf("%s/%s_SIC_contlimit_%s_%s_m%d_%s_cont.csv",        opt$plotfolder, opt$mode, opt$kernel, opt$integraltype, opt$bmass, errstring))
