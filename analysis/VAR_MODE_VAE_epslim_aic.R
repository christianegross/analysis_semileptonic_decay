library("hadron")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args)==1)
mode <- args[1]
stopifnot(mode=="DG" || mode=="DM" || mode=="DM2")

if(mode=="DG") {
    zlist <- c(0, 1, 2, 3)
}

if(mode=="DM") {
    zlist <- c(0, 1, 2, 3, 4)
}

if(mode=="DM2") {
    zlist <- c(0, 1, 2, 3, 4, 5)
}
fnfour <- function(par, x, boot.R, ...) par[1] + par[2] * x^2 + par[3] * x^4
fnsix <- function(par, x, boot.R, ...) par[1] + par[2] * x^2 + par[3] * x^4 + par[4] * x^6
epslist <- 1:20
par.guess <- rep(1, 3)

res <- data.frame(channel=c(), kernel=c(), theta=c(), w=c(), errtype=c(), iz=c(), DGDq2=c(), dDGDq2=c())

index <- 1

savefolder <- "tables_fnfour_20"
dividemass <- F
comment <- "_aic"
if(dividemass) comment <- "_dividemass"
pdf(sprintf("plots/%s_VAE_epslim%s.pdf", mode, comment), title="")

reslist <- list(channel=c(), kernel=c(), theta=c(), w=c(), errtype=c(), iz=c(), DGDq2=c(), bsDGDq2=array(NA, dim=c(1000, 2*2*10*4*length(zlist))))
for(channel in c("cd", "cs")) {
  for(kernel in c("sigmoid", "erf")) {
#~ for(channel in c("cd")) {
#~   for(kernel in c("erf")) {
    for(theta in as.character(c(1:9, 9.5))) {
      dat <- read.table(sprintf("%s/%s_VAE_contlim_%s_%s_th%s%s.csv", savefolder, mode, channel, kernel, theta, comment), header=TRUE)
      datlist <- readRDS(sprintf("%s/%s_VAE_contlim_%s_%s_th%s%s.RDS", savefolder, mode, channel, kernel, theta, comment))
      for(errtype in c("stat", "sys", "vol", "tot")) {
        for(iz in zlist) {
          dattmp <- dat[dat$errtype==errtype & dat$iz==iz, ]
          bs <- datlist$dat[, datlist$iz==iz & datlist$errtype==errtype]
          # print(dattmp)
          nafraction <- sum(is.na(bs)) / length(bs)
          
          myfit <- try(bootstrap.nlsfit(fn=fnfour, x=dattmp$eps, y=dattmp$DGDq2, bs=bs, par.guess=par.guess, na.rm=T, mask=epslist))
          if(!inherits(myfit, "try-error") & nafraction < 0.05) {
          plot(myfit, xlab="epsilon/mH", ylab=sprintf("%s/Dw^2", mode), plot.range=c(-0.1, max(myfit$x)*1.1), xlim=c(c(0, max(myfit$x))),
               main=paste("channel", channel, "kernel", kernel, "th", theta, "iz", iz, "err", errtype))
          plotwitherror(x=0, y=myfit$t0[1], dy=myfit$se[1], col=2, rep=T)
          newline <- data.frame(channel=channel, kernel=kernel, theta=theta, w=dattmp$w[1], errtype=errtype, 
                                iz=iz, DGDq2=myfit$t0[1], dDGDq2=myfit$se[1])
          res <- rbind(res, newline)
          reslist$channel[index]   <- channel
          reslist$kernel[index]    <- kernel          
          reslist$theta[index]     <- theta
          reslist$w[index]         <- dattmp$w[1]
          reslist$errtype[index]   <- errtype          
          reslist$iz[index]        <- iz
          reslist$DGDq2[index]     <- myfit$t0[1]          
          reslist$bsDGDq2[, index] <- myfit$t[, 1]
          } else {
              plotwitherror(x=dattmp$eps, y=dattmp$DGDq2, dy=apply(X=bs, MARGIN=2, FUN=sd), xlab="epsilon/mH", ylab=sprintf("%s/Dw^2", mode),
               main=paste("channel", channel, "kernel", kernel, "th", theta, "iz", iz, "err", errtype))
          newline <- data.frame(channel=channel, kernel=kernel, theta=theta, errtype=errtype, 
                                iz=iz, DGDq2=NA, dDGDq2=NA)
          res <- rbind(res, newline)
          reslist$channel[index]   <- channel
          reslist$kernel[index]    <- kernel          
          reslist$theta[index]     <- theta
          reslist$w[index]         <- NA
          reslist$errtype[index]   <- errtype          
          reslist$iz[index]        <- iz
          reslist$DGDq2[index]     <- NA          
          reslist$bsDGDq2[, index] <- rep(NA, 1000)
          }
          index <- index+1
        }
      }
    }
  }
}
print(res)
write.table(x=res, file=sprintf("%s/%s_VAE_epslim%s.csv", savefolder, mode, comment), col.names=T, row.names=F)
saveRDS(object=reslist, file=sprintf("%s/%s_VAE_epslim%s.RDS", savefolder, mode, comment))
