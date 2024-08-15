library("hadron")
fnfour <- function(par, x, boot.R, ...) par[1] + par[2] * x^2 + par[3] * x^4
zlist <- 0:4

res <- data.frame(channel=c(), kernel=c(), theta=c(), errtype=c(), iz=c(), DMDq2=c(), dDMDq2=c())

index <- 1

savefolder <- "tables"
dividemass <- T
comment <- ""
if(dividemass) comment <- "_dividemass"
pdf(sprintf("plots/DM_VAE_epslim%s.pdf", comment), title="")

reslist <- list(channel=c(), kernel=c(), theta=c(), errtype=c(), iz=c(), DMDq2=c(), bsDMDq2=array(NA, dim=c(1000, 2*2*10*4*4)))
for(channel in c("cd", "cs")) {
  for(kernel in c("sigmoid", "erf")) {
    for(theta in as.character(c(1:9, 9.5))) {
      dat <- read.table(sprintf("%s/DM_VAE_contlim_%s_%s_th%s%s.csv", savefolder, channel, kernel, theta, comment), header=TRUE)
      datlist <- readRDS(sprintf("%s/DM_VAE_contlim_%s_%s_th%s%s.RDS", savefolder, channel, kernel, theta, comment))
      for(errtype in c("stat", "sys", "vol", "tot")) {
        for(iz in zlist) {
          dattmp <- dat[dat$errtype==errtype & dat$iz==iz, ]
          bs <- datlist$dat[, datlist$iz==iz & datlist$errtype==errtype]
          # print(dattmp)
          
          myfit <- bootstrap.nlsfit(fn=fnfour, x=dattmp$eps, y=dattmp$DMDq2, bs=bs, par.guess=c(1, 1, 1))
          plot(myfit, xlab="epsilon/mH", ylab="DM/Dw^2", plot.range=c(-0.1, max(myfit$x)*1.1), xlim=c(c(0, max(myfit$x))),
               main=paste("channel", channel, "kernel", kernel, "th", theta, "iz", iz, "err", errtype))
          plotwitherror(x=0, y=myfit$t0[1], dy=myfit$se[1], col=2, rep=T)
          newline <- data.frame(channel=channel, kernel=kernel, theta=theta, errtype=errtype, 
                                iz=iz, DMDq2=myfit$t0[1], dDMDq2=myfit$se[1])
          res <- rbind(res, newline)
          reslist$channel[index]   <- channel
          reslist$kernel[index]    <- kernel          
          reslist$theta[index]     <- theta
          reslist$errtype[index]   <- errtype          
          reslist$iz[index]        <- iz
          reslist$DMDq2[index]     <- myfit$t0[1]          
          reslist$bsDMDq2[, index] <- myfit$t[, 1]

          
          
          print(newline)
          index <- index+1
        }
      }
    }
  }
}
print(res)
write.table(x=res, file=sprintf("%s/DM_VAE_epslim%s.csv", savefolder, comment), col.names=T, row.names=F)
saveRDS(object=reslist, file=sprintf("%s/DM_VAE_epslim%s.RDS", savefolder, comment))
