
source("/home/gross/Documents/heavymesons/scripts/extrapolate_sigma_templates/calc_DGDq2.R")
library("hadron")

erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

## extract the data for a single z, epsilon combo from the DGDq2 measurements
extractdata <- function(data, iz=0, icomb=0, iset=0, epsilons=seq(1, 19)-1) {
  namesel <- paste0("id", indexDG(inorm = 0, ieps = epsilons, icomb = icomb, idg = iz, nnorms = 1, neps = data[[2+iset]]$metadata$neps, NCOMBS = 5), 
                    "ieps", epsilons, "icomb", icomb, "iz", iz)
  # namesel <- paste0("id", ((iz)*5 + icomb)*data[[2+iset]]$metadata$neps + (epsilons), "ieps", epsilons, "icomb", icomb, "iz", iz)
  # print(namesel)
  # print(names(data[[2+iset]]))
  # print(names(data[[2+iset]][namesel]))
  # print(sapply(X=data[[2+iset]][namesel], FUN=getElement, name="DGDq2mean"))
  # print(sapply(X=data[[2+iset]][namesel], FUN=getElement, name="DGDq2mean")[1, ])
  DGamma <- unlist(sapply(X=data[[2+iset]][namesel], FUN=getElement, name="DGDq2mean")[1, ], use.names=F)
  dDGamma_sys <- unname(unlist(sapply(X=data[[2+iset]][namesel], FUN=getElement, name="sys")[1, ], use.names=F))
  dDGamma <- unlist(sapply(X=data[[2+iset]][namesel], FUN=getElement, name="DGDq2sd")[1, ], use.names=F)
  DGammaboot <- sapply(X=data[[2+iset]][namesel], FUN=getElement, name="DGDq2boots")[1, ]
  
  bs <- array(unlist(DGammaboot, use.names=F), dim=c(data[[2+iset]][[namesel[1]]]$bootnumber[[1]], length(epsilons)))
  x <- data[[2+iset]]$epsilons[epsilons]
  
  ## determine min A/A0_ref
  AA0ar <- array(NA, dim=c(data[[2+iset]]$metadata$nnorm, length(epsilons)))
  for (norm in seq(1, data[[2+iset]]$metadata$nnorm)) {
    A0ABCW_ref <- sapply(X=data[[2+iset]][namesel], FUN=getElement, name="A0ABCW_ref")[1, ]
    A0ABCW_ref_ar <- t(array(unlist(A0ABCW_ref, use.names=F), dim=c(6, length(epsilons))))
    AA0ar[norm, ] <- A0ABCW_ref_ar[, 2] / A0ABCW_ref_ar[, 1]
  }
  AA0 <- apply(X=AA0ar, MARGIN=2, FUN=min)
  # print(list(epsilons=x, DGammamean=DGamma, dDGamma=dDGamma, dDGamma_sys=dDGamma_sys, AA0=AA0))
  return(list(epsilons=x, DGammamean=DGamma, dDGamma=dDGamma, dDGamma_sys=dDGamma_sys, boot=bs, AA0=AA0))
  
}

# determine pull factor if one set of measurements is in the DGDq2.bin file format and the other has been interpolated and saved in a separate type of list
## https://arxiv.org/abs/2212.08467


determinesyserrfinitevolumeextrapolated <- function(resultpath1, filenames1, nerr1, L1, resultpath2, filenames2, nerr2, L2, th, 
                                                    savename, isets=c(0), icomb=c(0), mode="DG") {
  ## detect some errors
  stopifnot(length(filenames1) == length(filenames2))
  stopifnot(length(filenames1) == length(nerr1))
  stopifnot(length(filenames1) == length(nerr2))
  if (any(nerr1 != nerr2)) { print("WARNING: nerr are not the same, assume nerr1 everywhere")}
  if(mode=="DG") {
    NDG <- 4
    zlist=c(0, 1, 2, 3)
  }
  else if(mode=="DM")  {
    NDG <- 5
    zlist=c(0, 1, 2, 3, 4)
  }
  else stop("mode has to be DG or DM")
  
  
  result <- data.frame(w=NA, nerr=NA, iz=NA, icomb=NA, th = NA, epsilon=NA, P=NA, Delta=NA, L1=NA, L2=NA, dDGamma2=NA)
  
  # resultdat <- list(tsnk=c(), Nt=c(), w=c(), nerr=c(), iz=c(), icomb=c(), th=c(), includesys = c(), DGDq2 = c(), dDGDq2 = c(), dat = array(NA, dim=c(1000, length(tsnk)*2*length(isets)*length(icomb)*4)))
  parindex <- 1
  for (index in seq(1, length(filenames1))){
    data1 <- try(read_in_DGDq2(filename=filenames1[index], write=FALSE, resultpath=resultpath1[index], NDG = 5))
    print(paste0(resultpath2[index], "/", filenames2[index]))
    data2 <- try(readRDS(file=paste0(resultpath2[index], "/", filenames2[index])))
    if(!inherits(x=data1, what="try-error") && !inherits(x=data2, what="try-error")) {
      for (iset in isets) {
        mymasks <- list(c(NA))
        title <- paste("iset", iset, "th", th[index])
        print(title)
        for (iz in zlist) {
          for(icomb in c(0)) {
            
            
            unequallength <- FALSE
            minlength <- data1[[2+iset]]$metadata$neps
            # stopifnot(length(DGamma1) == length(DGamma2))
            if(data1[[2+iset]]$metadata$neps != data2$neps) {
              print("WARNING: DGamma1 and DGamma2 do not have the same length")
              unequallength <- TRUE
              maxlength <- max(data1[[2+iset]]$metadata$neps, data2$neps)
              minlength <- min(data1[[2+iset]]$metadata$neps, data2$neps)
            }
            
            ## determine which part of data to use, extract needed elements
            namesel1 <- paste0("id", ((iz)*5 + icomb)*data1[[2+iset]]$metadata$neps + (0:(minlength-1)), 
                               "ieps", (0:(minlength-1)), "icomb", icomb, "iz", iz)
            namesel2 <- paste0("icomb", icomb, "iz", iz, "epsilon", (0:(minlength-1)))
            
            
            
            ## read in and extract information from first ensemble
            
            DGamma1 <- unlist(sapply(X=data1[[2+iset]][namesel1], FUN=getElement, name="DGDq2mean")[1, ], use.names=F)
            dDGamma1 <- unlist(sapply(X=data1[[2+iset]][namesel1], FUN=getElement, name="DGDq2sd")[1, ], use.names=F)
            
            
            ## read in and extract information from second ensemble
            DGamma2 <- unlist(sapply(X=data2[["interpolate"]][namesel2], FUN=getElement, name="val"), use.names=F)
            dDGamma2 <- unlist(sapply(X=data2[["interpolate"]][namesel2], FUN=getElement, name="err"), use.names=F)
            
            names(DGamma2) <- NULL
            names(dDGamma2) <- NULL
            
            
            
            ## determine P according to eq. 40
            
            P <- (DGamma1 - DGamma2)/sqrt(dDGamma1^2+dDGamma2^2)
            
            
            ## determine Delta as in eq. 42
            
            Delta <- abs(DGamma1 - DGamma2) * erf(abs(P)/sqrt(2))
            # plot(x=data1[[2+iset]]$epsilons, y=Delta, xlab="epsilon", ylab="Delta", main=title)
            
            if(unequallength) {
              P[(minlength+1):maxlength] <- NA
              Delta[(minlength+1):maxlength] <- NA
              dDGamma2[(minlength+1):maxlength] <- NA
            }
            
            ## save result
            neps <- max(data1[[2+iset]]$metadata$neps, data2$neps)
            result <- rbind(result, data.frame(w=rep(data1[[2+iset]]$metadata$w, neps),
                                               nerr=rep(nerr1[index], neps), iz=rep(iz, neps), icomb=rep(icomb, neps),
                                               th=rep(th[index], neps), epsilon=data1[[2+iset]]$epsilons, P=P, Delta=Delta,
                                               L1=rep(L1[index], neps), L2=rep(L2[index], neps), dDGamma2=dDGamma2))
            
            
          }
        }
        plot(x=result$epsilon[result$th==th[index]][!is.na(result$P[result$th==th[index]])], y=result$P[result$th==th[index]][!is.na(result$P[result$th==th[index]])], 
             xlab="epsilon", ylab="P", main=title,
             col=result$iz[result$th==th[index]][!is.na(result$P[result$th==th[index]])]+1, type="p")
        lines(x=c(-1, 2), y=rep(0, 2))
        if (mode=="DG") legend(x="topleft", legend=c(0, 1, 2, "3   Z"), col=seq(1, 4), pch=c(1, 4), horiz=TRUE)
        if (mode=="DM") legend(x="topleft", legend=c(0, 1, 2, 3, "4   Z"), col=seq(1, 5), pch=c(1, 5), horiz=TRUE)
      }
    }
  }
  result <- result[-1, ]
  result <- result[ ! (result$iz==0 & (result$icomb==3 | result$icomb==4)), ]
  print(warnings())
  write.table(result, paste0(savename, ".csv"), row.names=F, col.names=T)
  return(invisible(result))
}

## determine pull factor if both datasets are in the DGDq2.bin format
## https://arxiv.org/abs/2212.08467

determinesyserrfinitevolume <- function(resultpath1, filenames1, nerr1, L1, resultpath2, filenames2, nerr2, L2, th, 
                                        savename, isets=c(0), icomb=c(0), mode="DG") {
  ## detect some errors
  stopifnot(length(filenames1) == length(filenames2))
  stopifnot(length(filenames1) == length(nerr1))
  stopifnot(length(filenames1) == length(nerr2))
  if(mode=="DG") {
    NDG <- 4
    zlist=c(0, 1, 2, 3)
  }
  else if(mode=="DM")  {
    NDG <- 5
    zlist=c(0, 1, 2, 3, 4)
  }
  else stop("mode has to be DG or DM")
  
  
  result <- data.frame(w=NA, nerr=NA, iz=NA, icomb=NA, th = NA, epsilon=NA, P=NA, Delta=NA, L1=NA, L2=NA, dDGamma2=NA)
  
  # resultdat <- list(tsnk=c(), Nt=c(), w=c(), nerr=c(), iz=c(), icomb=c(), th=c(), includesys = c(), DGDq2 = c(), dDGDq2 = c(), dat = array(NA, dim=c(1000, length(tsnk)*2*length(isets)*length(icomb)*4)))
  parindex <- 1
  for (index in seq(1, length(filenames1))){
    data1 <- try(read_in_DGDq2(filename=filenames1[index], write=FALSE, resultpath=resultpath1[index], NDG = NDG))
    data2 <- try(read_in_DGDq2(filename=filenames2[index], write=FALSE, resultpath=resultpath2[index], NDG = NDG))
    if(!inherits(x=data1, what="try-error") && !inherits(x=data2, what="try-error")) {
      for (iset in isets) {
        mymasks <- list(c(NA))
        title <- paste("iset", iset, "th", th[index])
        print(title)
        for (iz in zlist) {
          for(icomb in c(0)) {
            
            
            unequallength <- FALSE
            minlength <- data1[[2+iset]]$metadata$neps
            # stopifnot(length(DGamma1) == length(DGamma2))
            if(data1[[2+iset]]$metadata$neps != data2[[2+iset]]$metadata$neps) {
              print("WARNING: DGamma1 and DGamma2 do not have the same length")
              unequallength <- TRUE
              maxlength <- max(data1[[2+iset]]$metadata$neps, data2[[2+iset]]$metadata$neps)
              minlength <- min(data1[[2+iset]]$metadata$neps, data2[[2+iset]]$metadata$neps)
            }
            neps <- max(data1[[2+iset]]$metadata$neps, data2[[2+iset]]$metadata$neps)
            
            ## determine which part of data to use, extract needed elements
            namesel1 <- paste0("id", ((iz)*5 + icomb)*data1[[2+iset]]$metadata$neps + (0:(minlength-1)), 
                               "ieps", (0:(minlength-1)), "icomb", icomb, "iz", iz)
            namesel2 <- paste0("id", ((iz)*5 + icomb)*data2[[2+iset]]$metadata$neps + (0:(minlength-1)), 
                               "ieps", (0:(minlength-1)), "icomb", icomb, "iz", iz)
            
            
            
            ## read in and extract information from first ensemble
            
            DGamma1 <- unlist(sapply(X=data1[[2+iset]][namesel1], FUN=getElement, name="DGDq2mean")[1, ], use.names=F)
            dDGamma1 <- unlist(sapply(X=data1[[2+iset]][namesel1], FUN=getElement, name="DGDq2sd")[1, ], use.names=F)
            
            
            ## read in and extract information from second ensemble
            
            DGamma2 <- unlist(sapply(X=data2[[2+iset]][namesel2], FUN=getElement, name="DGDq2mean")[1, ], use.names=F)
            dDGamma2 <- unlist(sapply(X=data2[[2+iset]][namesel2], FUN=getElement, name="DGDq2sd")[1, ], use.names=F)
            
            
            ## determine P according to eq. 40
            
            P <- (DGamma1 - DGamma2)/sqrt(dDGamma1^2+dDGamma2^2)
            
            
            ## determine Delta as in eq. 42
            
            Delta <- abs(DGamma1 - DGamma2) * erf(abs(P)/sqrt(2))
            # plot(x=data1[[2+iset]]$epsilons, y=Delta, xlab="epsilon", ylab="Delta", main=title)
            
            if(unequallength) {
              P[(minlength+1):maxlength] <- NA
              Delta[(minlength+1):maxlength] <- NA
              dDGamma2[(minlength+1):maxlength] <- NA
            }
            
            ## save result
            result <- rbind(result, data.frame(w=rep(data1[[2+iset]]$metadata$w, neps),
                                               nerr=rep(nerr[index], neps), iz=rep(iz, neps), icomb=rep(icomb, neps),
                                               th=rep(th[index], neps), epsilon=unique(c(data1[[2+iset]]$epsilons, data2[[2+iset]]$epsilons)), P=P, Delta=Delta,
                                               L1=rep(L1[index], neps), L2=rep(L2[index], neps), dDGamma2=dDGamma2))
            
            
          }
        }
        plot(x=result$epsilon[result$th==th[index]][!is.na(result$P[result$th==th[index]])], y=result$P[result$th==th[index]][!is.na(result$P[result$th==th[index]])], 
             xlab="epsilon", ylab="P", main=title,
             col=result$iz[result$th==th[index]][!is.na(result$P[result$th==th[index]])]+1, type="p")
        lines(x=c(-1, 2), y=rep(0, 2))
        if (mode=="DG") legend(x="topleft", legend=c(0, 1, 2, "3   Z"), col=seq(1, 4), pch=c(1, 4), horiz=TRUE)
        if (mode=="DM") legend(x="topleft", legend=c(0, 1, 2, 3, "4   Z"), col=seq(1, 5), pch=c(1, 5), horiz=TRUE)
      }
    }
  }
  result <- result[-1, ]
  result <- result[ ! (result$iz==0 & (result$icomb==3 | result$icomb==4)), ]
  print(warnings())
  write.table(result, paste0(savename, ".csv"), row.names=F, col.names=T)
  return(invisible(result))
}

