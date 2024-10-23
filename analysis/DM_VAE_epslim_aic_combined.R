library("hadron")
source("/hiskp4/gross/heavymesons/helpscripts/combindedfit.R")
fnfour <- function(par, x, boot.R, ...) par[1] + par[2] * x^2 + par[3] * x^4
fnsix <- function(par, x, boot.R, ...) par[1] + par[2] * x^2 + par[3] * x^4 + par[4] * x^6
zlist <- 0:4

fnfourcomb <- function (par, x, boot.r, maskfn, ...) {
  return (par[1]  +  par[2] * x^2 * maskfn + par[3] * x^4 * maskfn + par[4] * x^2 * (!maskfn) + par[5] * x^4 * (!maskfn))
}

res <- data.frame(channel=c(), kernel=c(), theta=c(), errtype=c(), iz=c(), DGDq2=c(), dDGDq2=c())

index <- 1

savefolder <- "tables_fnfour_12"
dividemass <- F
comment <- "_aic"
if(dividemass) comment <- "_dividemass"
pdf(sprintf("plots/DM_VAE_epslim%s_combined.pdf", comment), title="")

reslist <- list(channel=c(), kernel=c(), theta=c(), errtype=c(), iz=c(), DGDq2=c(), bsDGDq2=array(NA, dim=c(1000, 2*2*10*4*length(zlist))))
for(channel in c("cd", "cs")) {
#~ for(channel in c("cd")) {
#~   for(kernel in c("erf")) {
    for(theta in as.character(c(1:9, 9.5))) {
    datsigmoid <- read.table(sprintf("%s/DM_VAE_contlim_%s_sigmoid_th%s%s.csv", savefolder, channel, theta, comment), header=TRUE)
    datlistsigmoid <- readRDS(sprintf("%s/DM_VAE_contlim_%s_sigmoid_th%s%s.RDS", savefolder, channel, theta, comment))
    daterf <- read.table(sprintf("%s/DM_VAE_contlim_%s_erf_th%s%s.csv", savefolder, channel, theta, comment), header=TRUE)
    datlisterf <- readRDS(sprintf("%s/DM_VAE_contlim_%s_erf_th%s%s.RDS", savefolder, channel, theta, comment))
      for(errtype in c("stat", "sys", "vol", "tot")) {
        for(iz in zlist) {
        dattmpsigmoid <- datsigmoid[datsigmoid$errtype==errtype & datsigmoid$iz==iz, ]
        bssigmoid <- datlistsigmoid$dat[, datlistsigmoid$iz==iz & datlistsigmoid$errtype==errtype]
        dattmperf <- daterf[daterf$errtype==errtype & daterf$iz==iz, ]
        bserf <- datlisterf$dat[, datlisterf$iz==iz & datlisterf$errtype==errtype]
        # print(dattmp)a
        
        bscomb <- array(c(bssigmoid, bserf), dim=c(length(bssigmoid[, 1]), length(bssigmoid[1, ]) + length(bserf[1, ])))
        
        maskfit <- c(rep(T, 12), rep(F, length(dattmpsigmoid$eps)-12), rep(T, 12), rep(F, length(dattmperf$eps)-12))
        maskfn <- c(rep(T, length(dattmpsigmoid$eps)), rep(F, length(dattmperf$eps)))
        myfit <- try(bootstrap.nlsfit(fn=fnfourcomb, x=c(dattmpsigmoid$eps, dattmperf$eps), 
                                      y=c(dattmpsigmoid$DGDq2, dattmperf$DGDq2), bs=bscomb, 
                                      par.guess=c(1, 1, 1,1, 1), mask=maskfit, 
                                      maskfn=maskfn[maskfit], na.rm=T))
        myfit$maskfn <- maskfn
        if(!inherits(myfit, "try-error")) {
          plotcombined(myfit, xlab="epsilon/mH", ylab="DM/Dw^2", xlim=c(c(0, max(myfit$x))),
               main=paste("channel", channel, "th", theta, "iz", iz, "err", errtype))
          plotwitherror(x=0, y=myfit$t0[1], dy=myfit$se[1], col=2, rep=T)
          newline <- data.frame(channel=channel, kernel="combined", theta=theta, errtype=errtype, 
                                iz=iz, DGDq2=myfit$t0[1], dDGDq2=myfit$se[1])
          res <- rbind(res, newline)
          reslist$channel[index]   <- channel
          reslist$kernel[index]    <- "combined"          
          reslist$theta[index]     <- theta
          reslist$errtype[index]   <- errtype          
          reslist$iz[index]        <- iz
          reslist$DGDq2[index]     <- myfit$t0[1]          
          reslist$bsDGDq2[, index] <- myfit$t[, 1]
          }

          
          
          print(newline)
          index <- index+1
        }
      }
    }
  }

print(res)
write.table(x=res, file=sprintf("%s/DM_VAE_epslim%s_combined.csv", savefolder, comment), col.names=T, row.names=F)
saveRDS(object=reslist, file=sprintf("%s/DM_VAE_epslim%s_combined.RDS", savefolder, comment))
