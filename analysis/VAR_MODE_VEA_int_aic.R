source("/hiskp4/gross/heavymesons/helpscripts/integrate.R")
source("/hiskp4/gross/heavymesons/helpscripts/splineintegration_functions.R")
library("hadron")
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

errlist <- c("stat", "sys", "vol", "tot")
doplot <- T
fnlin <- function(par, x, boot.R, ...) par[1] + par[2] * x
fncon <- function(par, x, boot.R, ...) par[1] + 0*x


savefolder <- "tables_fnfour_20"


upperboundarycd <- 0.8724
upperboundarycs <- 0.7770
channels <- c("cd", "cs")
channelboundaries <- c(upperboundarycd, upperboundarycs)
continue <- c(T, F)
replacelower <- c(F, T)
replaceindex <- c(0, 10)
confcounter <- 1

pdf(sprintf("plots/%s_VEA_int_%s.pdf", mode, comment), title="")

for(channel_index in seq_along(channels)) {
  channel <- channels[channel_index]
  upperbound <- channelboundaries[channel_index]
    for(kernel in c("sigmoid", "erf")) {

contlim <- readRDS(sprintf("%s/%s_VEA_contlim_%s_%s_%s.RDS", savefolder, mode, channel, kernel, comment))

contlimtable <- read.table(sprintf("%s/%s_VEA_contlim_%s_%s_%s.csv", savefolder, mode, channel, kernel, comment))



enslist <- list(contlim)
tablelist <- list(contlimtable)
names <- c("contlim")


res <- data.frame(name=NA, ensno=NA, errtype=NA, iz=NA, intspline=NA, dintspline=NA, intbsspline=NA, inttrap=NA, dinttrap=NA, intbstrap=NA)
reslist <- list()


confcounter <- 1
thetas <- c(seq(1, 9), 9.5)
for (i in seq_along(enslist)){
  ens <- enslist[[i]]
  mytable <- tablelist[[i]]
  for (errtype in errlist) {
    for (iz in zlist) {
      title <- sprintf("ens %s iz %d errtype %s", names[i], iz, errtype)
      titlelist <- sprintf("ens%s-iz%d-errtype%s", names[i], iz, errtype)
      print(title)
      ## we include the error on a in the continuum limit
      bsamples <- array(rep(0, 11000), dim=c(1000, 11))
      y <- c(0)
      dy <- c(0)
      x <- c(0)
      dx <- c(0)
      for (index in seq_along(thetas)) {
        th <- thetas[index]
        if(!is.na(mytable$q[abs(mytable$th - th) < 1e-2][1]^2) && names[i] != "contlim") {
          mask <- abs(ens$th - th) < 1e-2 & ens$iz == iz & ens$errtype==errtype
          mask[is.na(mask)] <- F
          masktable <- abs(mytable$th - th) < 1e-2 & mytable$iz == iz & mytable$errtype==errtype
          masktable[is.na(masktable)] <- F
          y[index+1] <- ens$DGDq2gev[mask]
          dy[index+1] <- ens$dDGDq2gev[mask]
          try(bsamples[, index+1] <- ens$datgev[, mask])
          x[index+1] <- mytable$q[masktable]^2
        } else if(!is.na(mytable$q[abs(mytable$th - th) < 1e-2][1]^2) && names[i] == "contlim") {
          mask <- abs(ens$th - th) < 1e-2 & ens$iz == iz & ens$errtype==errtype
          mask[is.na(mask)] <- F
          masktable <- abs(mytable$th - th) < 1e-2 & mytable$iz == iz & mytable$errtype==errtype
          masktable[is.na(masktable)] <- F
          y[index+1] <- ens$fit[[which(mask)]]$weighted$mean
          dy[index+1] <- ens$fit[[which(mask)]]$weighted$sd
          try(bsamples[, index+1] <- ens$fit[[which(mask)]]$weighted$bootsamples)
          x[index+1] <- mytable$q[masktable]^2
        } else {
          print(sprintf("failure for th %s", index))
        }
      }
      
      spline <- interpSpline(x, y)
      len <- length(x)
      slopebs <- (bsamples[, len] - bsamples[, len-1])/(x[len] - x[len-1])
      yupperbs <- bsamples[, len] + (upperboundarycd-x[len]) * slopebs
      
      plotwitherror(x=x, y=y, dy=dy, xlab="q^2", ylab=sprintf("%sDq^2", mode), main=paste(channel, kernel, title))
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
      res <- rbind(res, data.frame(name=names[i], ensno=i, errtype=errtype, iz=iz, 
                                   intspline=meanintspline, dintspline=sd(bsintspline, na.rm=T), intbsspline=mean(bsintspline, na.rm=T), 
                                   inttrap=meaninttrap, dinttrap=sd(bsinttrap, na.rm=T), intbstrap=mean(bsinttrap, na.rm=T)))
      
      reslist[[paste0(names[i], errtype, "iz", iz, "bsintspline")]] <- bsintspline
      reslist[[paste0(names[i], errtype, "iz", iz, "spline")]] <- spline
      reslist[[paste0(names[i], errtype, "iz", iz, "bsinttrap")]] <- bsinttrap
      confcounter <- confcounter + 1
    }
  }
}

res <- res[-1, ]
res
write.table(x=res, file=sprintf("%s/%s_VEA_%s_%s_int_%s.csv", savefolder, mode, channel, kernel, comment), row.names=F, col.names=T)

reslist$info <- res
saveRDS(object=reslist, file=sprintf("%s/%s_VEA_%s_%s_int_%s.RDS", savefolder, mode, channel, kernel, comment)) 

for(errtype in errlist) {
  mask <- res$errtype==errtype
  plotwitherror(x=res$iz[mask]+ 0.1*(res$ensno[mask]-1), y=res$intspline[mask], dy=res$dintspline[mask], col=res$ensno[mask], pch=res$ensno[mask], 
                main=paste("err =", errtype), xlab="Z", ylab=sprintf("%s [a.u.]", mode))
  legend(legend=names, x="topleft", col=seq(1, max(res$ensno)), pch=seq(1, max(res$ensno)))
}

}
}

dev.off()
