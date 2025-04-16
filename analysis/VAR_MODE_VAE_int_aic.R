library("hadron")
library("splines")
comment <- ""
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 2) {
    comment <- args[1]
    mode <- args[2]
}
if(length(args) != 2) stop("pass the comment for the calculations and the mode")
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

kernels <- c("sigmoid", "erf")
if( grepl("combined", comment, fixed=TRUE)) kernels <- c("combined")


source("/hiskp4/gross/heavymesons/helpscripts/integrate.R")
source("/hiskp4/gross/heavymesons/helpscripts/splineintegration_functions.R")

errlist <- c("stat", "sys", "vol", "tot")
dividemass <- F
if(dividemass) comment <- "_dividemass"
savefolder <- "tables_fnfour_20"

mytable <- read.table(sprintf("%s/%s_VAE_epslim%s.csv", savefolder, mode, comment), header=TRUE)
mydata <- readRDS(sprintf("%s/%s_VAE_epslim%s.RDS", savefolder, mode, comment))
print(sprintf("%s/%s_VAE_epslim%s.RDS", savefolder, mode, comment))
res <- data.frame(channel=NA, kernel=NA, ensno=NA, errtype=NA, iz=NA, intspline=NA, dintspline=NA, intbsspline=NA, inttrap=NA, dinttrap=NA, intbstrap=NA)
reslist <- list()
upperboundarycd <- 0.8724
upperboundarycs <- 0.7770


channels <- c("cd", "cs")
channelboundaries <- c(upperboundarycd, upperboundarycs)
continue <- c(T, F)
replacelower <- c(F, T)
replaceindex <- c(0, 10)
confcounter <- 1

pdf(sprintf("plots/%s_VAE_int%s.pdf", mode, comment), title="")

m_Ds <- 1.96835
print(zlist)
for(channel_index in seq_along(channels)) {
  channel <- channels[channel_index]
  upperbound <- channelboundaries[channel_index]
  for(kernel in kernels) {
    for (errtype in errlist) {
      for (iz in zlist) {
        title <- sprintf("iz %d errtype %s", iz, errtype)
        titlelist <- sprintf("iz%d-errtype%s", iz, errtype)
        print(title)
        ## we include the error on a in the continuum limit
        bsamples <- array(rep(0, 11000), dim=c(1000, 11))
        y <- c(0, mytable$DGDq2[mytable$iz==iz & mytable$errtype==errtype & mytable$channel==channel & mytable$kernel==kernel])
        dy <- c(0, mytable$dDGDq2[mytable$iz==iz & mytable$errtype==errtype & mytable$channel==channel & mytable$kernel==kernel])
        x <- c(0, mytable$w[mytable$iz==iz & mytable$errtype==errtype & mytable$channel==channel & mytable$kernel==kernel])^2*m_Ds^2
        bsamples[, 2:11] <- mydata$bsDGDq2[, which(mydata$iz==iz & mydata$errtype==errtype & mydata$channel==channel & mydata$kernel==kernel)]
        
        print(x)
        print(y)
        
        spline <- interpSpline(x, y)
        len <- length(x)
        slopebs <- (bsamples[, len] - bsamples[, len-1])/(x[len] - x[len-1])
        yupperbs <- bsamples[, len] + (upperboundarycd-x[len]) * slopebs
        
        plotwitherror(x=x, y=y, dy=dy, xlab="w^2", ylab=sprintf("%s/Dw^2", mode), main=paste(channel, kernel, title))
        xval <- seq(min(x), upperboundarycd, length.out=500)
        lines(x=xval, y=predict(object=spline, x=xval)$y, col="red", lty=2)
        lines(x, y, col="blue", lty=3)
        # plotwitherror(x=upperboundarycd, y=mean(yupperbs), dy=sd(yupperbs), col="blue", rep=T)
        lines(x=c(x[len], upperboundarycd), y=c(y[len], mean(yupperbs)), col="blue", lty=3)
        plotwitherror(x=x, y=y, dy=dy, rep=TRUE)
        legend(x="topright", legend=c("meas", "spline", "trapezoidal"), col=c(1, "red", "blue"), pch=c(1, NA, NA), lty=c(NA, 2, 3))
        
        meanintspline <- splineintegral(yval=y, xval=x, continue = continue[channel_index], 
                                        replacelower=replacelower[channel_index], higherlimit = upperbound, 
                                        lowerlimit=upperbound, replaceindex=replaceindex[channel_index])
        meaninttrap <- trapezoidal(yval=y, xval=x, continue = continue[channel_index], 
                                   replacelower=replacelower[channel_index], higherlimit = upperbound, 
                                   lowerlimit=upperbound, replaceindex=replaceindex[channel_index])
        bsintspline <- apply(X=bsamples, MARGIN=1, FUN=splineintegral, xval=x, continue = continue[channel_index], 
                             replacelower=replacelower[channel_index], higherlimit = upperbound, 
                             lowerlimit=upperbound, replaceindex=replaceindex[channel_index])
        bsinttrap <- apply(X=bsamples, MARGIN=1, FUN=trapezoidal, xval=x, continue = continue[channel_index], 
                           replacelower=replacelower[channel_index], higherlimit = upperbound, 
                           lowerlimit=upperbound, replaceindex=replaceindex[channel_index])
        res <- rbind(res, data.frame(channel=channel, kernel=kernel, ensno=confcounter, errtype=errtype, iz=iz, 
                                     intspline=meanintspline, dintspline=sd(bsintspline, na.rm=T), intbsspline=mean(bsintspline, na.rm=T), 
                                     inttrap=meaninttrap, dinttrap=sd(bsinttrap, na.rm=T), intbstrap=mean(bsinttrap, na.rm=T)))
        
        reslist[[paste0(channel, kernel, errtype, "iz", iz, "bsintspline")]] <- bsintspline
        reslist[[paste0(channel, kernel, errtype, "iz", iz, "spline")]] <- spline
        reslist[[paste0(channel, kernel, errtype, "iz", iz, "bsinttrap")]] <- bsinttrap
        confcounter <- confcounter + 1
      }
    }
  }
}
res <- res[-1, ]
res
write.table(x=res, file=sprintf("%s/%s_VAE_int%s.csv", savefolder, mode, comment), row.names=F, col.names=T)

reslist$info <- res
saveRDS(object=reslist, file=sprintf("%s/%s_VAE_int%s.RDS", savefolder, mode, comment)) 

dev.off()
