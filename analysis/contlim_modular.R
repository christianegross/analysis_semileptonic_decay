library("hadron")

library("optparse")


if (TRUE) {
  # set option list
  option_list <- list(
    make_option(c("-t", "--momentum"), type = "character", default = "th2",
                help = "momentum [default %default]"),
    make_option(c("-m", "--mode"), type = "character", default = "DG",
                help = "mode: DG, DM, DM2 [default %default]"),
    make_option(c("-k", "--kernel"), type = "character", default = "sigmoid",
                help = "kernel [default %default]"),
    make_option(c("-f", "--folder"), type = "character", default = "-1",
                help = "folder with results [default %default]"),
    make_option(c("-s", "--subfolder"), type = "character", default = "-1",
                help = "subfolder after ensemble with results [default %default]"),
    make_option(c("-p", "--plotfolder"), type = "character", default = "-1",
                help = "folder where plots and tables are stored [default %default]"),
    make_option(c("-b", "--bmass"), type = "integer", default = "-1",
                help = "index for bmass [default %default]"),
    make_option(c("-e", "--error"), type = "character", default = "stat",
                help = "index for bmass [default %default]")
    
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


B64 <- readRDS(sprintf("%s/cB211.07.64/%s/%s/chosen_%sBJnu_m%d_%s_ik0_boots_err_%s.RDS", opt$folder, opt$momentum, opt$subfolder, ifelse(opt$kernel=="erf", "erf_", ""), opt$bmass, opt$momentum, opt$error))
nsigma <- dim(B64)[2]
nboot <- dim(B64)[3]

C80 <- readRDS(sprintf("%s/cC211.06.80/%s/%s/chosen_%sBJnu_m%d_%s_ik0_boots_err_%s.RDS", opt$folder, opt$momentum, opt$subfolder, ifelse(opt$kernel=="erf", "erf_", ""), opt$bmass, opt$momentum, opt$error))
stopifnot(nsigma == dim(C80)[2])
stopifnot(nboot == dim(C80)[3])

D96 <- readRDS(sprintf("%s/cD211.054.96/%s/%s/chosen_%sBJnu_m%d_%s_ik0_boots_err_%s.RDS", opt$folder, opt$momentum, opt$subfolder, ifelse(opt$kernel=="erf", "erf_", ""), opt$bmass, opt$momentum, opt$error))
stopifnot(nsigma == dim(D96)[2])
stopifnot(nboot == dim(D96)[3])


reslist <- list()

iz <- 1
isigma <- 1

fitfn1 <- function(pars, x, boot.r, ...) pars[1] + x*pars[2]
fitfn2 <- function(pars, x, boot.r, ...) pars[1] + x*0

pdf(sprintf("%s/contlimit_%s_m%d_%s_%s.pdf", opt$plotfolder, opt$kernel, opt$bmass, opt$momentum, opt$error))

for(iz in 0:zmax) {
  for(isigma in 0:(nsigma-1)) {
#~   for(isigma in 0:2) {
    boots <- array(c(B64[iz+1, isigma+1, ], C80[iz+1, isigma+1, ], D96[iz+1, isigma+1, ]), dim=c(nboot, 3))
    
    afm <- c(0.07957, 0.06821, 0.05692)
    fit1 <- bootstrap.nlsfit(x=afm^2, y=apply(boots, 2, mean), bsamples=boots, fn=fitfn1, par.guess=c(1, 1))
    weights <- c(exp(-1*(fit1$chisqr + 2 * 2 - length(fit1$x))/2))
    fit2 <- try(bootstrap.nlsfit(x=afm^2, y=apply(boots, 2, mean), bsamples=boots, fn=fitfn2, par.guess=c(1), mask=c(F, T, T))
    if(inherits(fit2, "try-error")){
        weights[2] <- 0
        fit2 <- list(t0=c(0), se=c(0), t=array(rep(0, nboot), dim=c(nboot, 1)), failed=T)
    } else {
        weights[2] <- c(exp(-1*(fit2$chisqr + 2 * 1 - length(fit1$x))/2))
        fit2$failed=F
    }
    weights <- weights/sum(weights)
    print(weights)
    average <- c(fit1$t0[1], fit2$t0[2])*weights
    averageboots <- fit1$t[, 1]*weights[1] + fit2$t[, 1]*weights[2]
    averagese <- sd(averageboots)
    
    layout.matrix <- matrix(c(1, 2, 2, 2), nrow = 1, ncol = 4)
    layout(mat = layout.matrix,
           heights = 1, # Heights of the two rows
           widths = 1) # Widths of the two columns

    mai <- par("mai")
    par(mai=c(mai[1], mai[2], mai[3], 0))
    
    myhist <- hist(averageboots, plot = FALSE, breaks=50)
    nbars <- length(myhist$counts)
    midbreaks <- (myhist$breaks[2:(nbars+1)] + myhist$breaks[1:nbars]) / 2
    plot(x=myhist$counts, y=midbreaks, 
         xlim=c(max(myhist$counts), min(myhist$counts)),type="l",
         ylim=range(myhist$breaks, fit1$t0[1] + c(-1, 1)*fit1$se[1], fit1$y+fit1$dy, fit1$y-fit1$dy),
         xaxt="n", xlab="", ylab="DGDq2")

    par(mai=c(mai[1], 0, mai[3], mai[2]))
    plot(fit1, 
         xlab="a^2[fm^2]", ylab="", main=paste("Z", iz, "sigma", isigma), 
         plot.range=c(0, 0.07957^2), xlim=c(0, 0.07957^2), 
         ylim=range(myhist$breaks, fit1$t0[1] + c(-1, 1)*fit1$se[1], fit1$y+fit1$dy, fit1$y-fit1$dy),
         yaxt="n", col="black", col.band="black", opacity.band=weights[1])
    if(!fit2$failed) {
        plot(fit2, rep=T,
         plot.range=c(0, 0.06821^2),
         yaxt="n", col="blue", col.band="blue", opacity.band=weights[2])
     }
    axis(side=4)
    par(mai=mai)
    
    pull <- (fit1$t0[1]-fit1$y[3])/fit1$se[1]
    dsys <- abs(fit1$t0[1]-fit1$y[3])*erf(abs(pull)/sqrt(2))
    
    res <- list(mean=fit1$t0[1], sd=fit1$se[1], boot=fit1$t[, 1], pull=pull, cutoff=(fit1$t0[1]-fit1$y[3])/fit1$t0[1], 
                dsys=dsys, bootsys=fit1$t[, 1] + rnorm(nboot, 0, dsys))
    reslist[[paste0("iz", iz, "sigma", isigma)]] <- res
  }
}

saveRDS(object=reslist, file=sprintf("%s/contlimit_%s_m%d_%s_%s.RDS", opt$plotfolder, opt$kernel, opt$bmass, opt$momentum, opt$error))
#~ reslist[[1]]
#~ unlist(sapply(reslist, getElement, name="pull"))
#~ reslist$metadata <- opt
#~ opt
