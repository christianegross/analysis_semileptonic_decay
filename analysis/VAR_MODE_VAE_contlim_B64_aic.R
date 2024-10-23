source("/hiskp4/gross/heavymesons/helpscripts/contlimitbeforeepsilonlimit_aic.R")
library("hadron")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args)==1)
mode <- args[1]
stopifnot(mode=="DG" || mode=="DM" || mode=="DM2")

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

errlist <- c("stat", "sys", "vol", "tot")
doplot <- T
fnlin <- function(par, x, boot.R, ...) par[1] + par[2] * x
fncon <- function(par, x, boot.R, ...) par[1] + 0*x

nerr <- rep(2, 1)

maxAmin <- 1

dividemass <- F
comment <- ""
if(dividemass) comment <- "_dividemass"
savefolder <- "tables_fnfour_20"
comment <- "_aic"
reslist <- list()
names <- c()
index <- 1


for(channel in c("cd", "cs")) {
    for(kernel in c("sigmoid", "erf")) {
for(theta in as.character(c(1:9, 9.5))) {
  pdf(sprintf("plots/%s_%s_%s_VAE_contlim_th%s%s.pdf", mode, channel, kernel, theta, comment),)
  savename <- sprintf("%s/%s_VAE_contlim_%s_%s_th%s%s", savefolder, mode, channel, kernel, theta, comment)
  files <- sprintf("%s_new/output%s/%s.bin", kernel, output, output)
  mine <- determineDGDq2_contlim_aic(resultpathlist=list(sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.64/th%s/", channel, theta),
                                                     sprintf("/hiskp4/gross/heavymesons/data/%s/cC211.06.80_600/th%s/", channel, theta),
                                                     sprintf("/hiskp4/gross/heavymesons/data/%s/cD211.054.96/th%s/", channel, theta),
                                                     sprintf("/hiskp4/gross/heavymesons/data/%s/cE211.044.112_300/th%s/", channel, theta)),
                                 filenames = c(files, files, files, files), 
                                 th = theta, nerr = nerr, amin = maxAmin, savename = savename,
                                 fitfnlist=list(fnlin, fncon), 
                                 par.guess=list(c(1, 1), c(1)), 
                                 legendargs=list(x="top", legend=c("linear", "constant"), ncol=2),
                                 errors=c("stat", "sys", "vol", "tot"), doplot=TRUE, 
                                 volumetablelist = c(sprintf("tables/volume_interpolations/%s_tryB64factor_%s_%s.csv", mode, channel, kernel), 
                                                     sprintf("tables/volume_interpolations/%s_tryB64factor_%s_%s.csv", mode, channel, kernel), 
                                                     sprintf("tables/volume_interpolations/%s_tryB64factor_%s_%s.csv", mode, channel, kernel), 
                                                     sprintf("tables/volume_interpolations/%s_tryB64factor_%s_%s.csv", mode, channel, kernel)),
                                 numberspacings = 4, neps=20, afm=c(0.07957, 0.06821, 0.05692, 0.04891), amds = c(0.8, 0.684, 0.57, 0.49), 
                                 dividemass = dividemass, zlist=zlist, maxz=maxz, NDG=NDG)
  dev.off()
  reslist[[index]] <- mine
  names[index] <- paste0("channel", channel, "kernel", kernel, "theta", theta)
  index <- index + 1
}
}
}
names(reslist) <- names
saveRDS(reslist, sprintf("%s/%s_VAE_contlim_fits%s.RDS", savefolder, mode, comment))

