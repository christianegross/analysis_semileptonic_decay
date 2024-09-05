
source("/hiskp4/gross/heavymesons/helpscripts/contlimitbeforeepsilonlimit.R")
library("hadron")
zlist <- c(4, 0, 1, 2, 3)
errlist <- c("stat", "sys", "vol", "tot")
doplot <- T
fnlin <- function(par, x, boot.R, ...) par[1] + par[2] * x
fncon <- function(par, x, boot.R, ...) par[1]

nerr <- rep(2, 1)

maxAmin <- 1

dividemass <- F
comment <- ""
if(dividemass) comment <- "_dividemass"
savefolder <- "tables_fnfour_12"

for(channel in c("cd", "cs")) {
    for(kernel in c("sigmoid", "erf")) {
for(theta in as.character(c(1:9, 9.5))) {
  pdf(sprintf("plots/DM_%s_%s_VAE_contlim_th%s%s.pdf", channel, kernel, theta, comment),)
  savename <- sprintf("%s/DM_VAE_contlim_%s_%s_th%s", savefolder, channel, kernel, theta)
  files <- sprintf("%s_new/outputDMDq2/DMDq2.bin", kernel)
  mine <- determineDGDq2_contlim(resultpathlist=list(sprintf("/hiskp4/gross/heavymesons/data/%s/cB211.07.64/th%s/", channel, theta),
                                                     sprintf("/hiskp4/gross/heavymesons/data/%s/cC211.06.80_600/th%s/", channel, theta),
                                                     sprintf("/hiskp4/gross/heavymesons/data/%s/cD211.054.96/th%s/", channel, theta),
                                                     sprintf("/hiskp4/gross/heavymesons/data/%s/cE211.044.112_300/th%s/", channel, theta)),
                                 filenames = c(files, files, files, files), 
                                 th = theta, nerr = nerr, amin = maxAmin, savename = savename,
                                 fitfn=fnlin, 
                                 par.guess=c(1, 1), 
                                 errors=c("stat", "sys", "vol", "tot"), doplot=TRUE, 
                                 volumetablelist = c(sprintf("tables/volume_interpolations/DM_tryB64factor_%s_%s.csv", channel, kernel), 
                                                     sprintf("tables/volume_interpolations/DM_tryB64factor_%s_%s.csv", channel, kernel), 
                                                     sprintf("tables/volume_interpolations/DM_tryB64factor_%s_%s.csv", channel, kernel), 
                                                     sprintf("tables/volume_interpolations/DM_tryB64factor_%s_%s.csv", channel, kernel)),
                                 numberspacings = 4, neps=15, afm=c(0.07957, 0.06821, 0.05692, 0.04891), amds = c(0.8, 0.684, 0.57, 0.49), 
                                 dividemass = dividemass, zlist=c(0, 1, 2, 3, 4), maxz=5, NDG=5)
  dev.off()
}
}
}
