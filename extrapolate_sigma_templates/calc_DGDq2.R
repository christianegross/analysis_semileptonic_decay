#~ library("hadron")
source("/home/gross/Documents/heavymesons/scripts_christiane/extrapolate_sigma_templates/read_in_binary.R")

determineDGDq2 <- function(resultpath, filenames, tsnk, Nt, th, 
        nerr, amin, savename, fitfn, par.guess, isets=c(0), icomb=c(0), epsuplim=0, epslowlim=0) {
    ## detect some errors
    stopifnot(length(filenames) == length(tsnk))
    stopifnot(length(filenames) == length(Nt))
    stopifnot(length(filenames) == length(th))
    stopifnot(length(filenames) == length(nerr))
    
    
    result <- data.frame(tsnk=NA, Nt=NA, w=NA, nerr=NA, iz=NA, icomb=NA, th = NA, DGDq2=NA, dDGDq2=NA, chi=NA, p=NA, includesys=NA, fitcorr=NA)
    
    resultdat <- list(tsnk=c(), Nt=c(), w=c(), nerr=c(), iz=c(), icomb=c(), th=c(), includesys = c(), DGDq2 = c(), dDGDq2 = c(), dat = array(NA, dim=c(1000, length(tsnk)*2*length(isets)*length(icomb)*4)))
    parindex <- 1
for (index in seq(1, length(filenames))){
    data <- read_in_DGDq2(filename=filenames[index], write=FALSE, resultpath=resultpath)
for (iset in isets) {
    mymasks <- list(c(NA))
    for (iz in c(0, 1, 2, 3)) {
        for(icomb in c(0)) {

            ## determine which part of data to use, extract needed elements
            namesel <- paste0("id", ((iz)*5 + icomb)*data[[2+iset]]$metadata$neps + (epslowlim:(data[[2+iset]]$metadata$neps-1-epsuplim)), "ieps", epslowlim:(data[[2+iset]]$metadata$neps-1-epsuplim), "icomb", icomb, "iz", iz)
            title=paste("iset", iset, "iz", iz, "icomb", icomb, "tsnk", tsnk[index], "Nt", Nt[index], "th", th[index], "nerr", nerr[index])
            print(title)


            DGamma <- unlist(sapply(X=data[[2+iset]][namesel], FUN=getElement, name="DGDq2mean")[1, ], use.names=F)
            dDGamma_sys <- unlist(sapply(X=data[[2+iset]][namesel], FUN=getElement, name="sys")[1, ], use.names=F)
            dDGamma <- sqrt(unlist(sapply(X=data[[2+iset]][namesel], FUN=getElement, name="DGDq2sd")[1, ], use.names=F))
            DGammaboot <- sapply(X=data[[2+iset]][namesel], FUN=getElement, name="DGDq2boots")[1, ]

            bs <- array(unlist(DGammaboot, use.names=F), dim=c(data[[2+iset]][[namesel[1]]]$bootnumber[[1]], data[[2+iset]]$metadata$neps-epslowlim-epsuplim))

            ## determine min A/A0_ref
            AA0ar <- array(rep(NA, data[[2+iset]]$metadata$neps*data[[2+iset]]$metadata$nnorm), dim=c(data[[2+iset]]$metadata$nnorm, data[[2+iset]]$metadata$neps-epslowlim-epsuplim))
            for (norm in seq(1, data[[2+iset]]$metadata$nnorm)) {
                A0ABCW_ref <- sapply(X=data[[2+iset]][namesel], FUN=getElement, name="A0ABCW_ref")[1, ]
                A0ABCW_ref_ar <- t(array(unlist(A0ABCW_ref, use.names=F), dim=c(6, data[[2+iset]]$metadata$neps-epslowlim-epsuplim)))
                AA0ar[norm, ] <- A0ABCW_ref_ar[, 2] / A0ABCW_ref_ar[, 1]
            }
            AA0 <- apply(X=AA0ar, MARGIN=2, FUN=min)

            ## determine masks for fit
            if(iz != 3) mymasks[[iz+1]] <- (AA0 < amin)
            if (iz==3) {
                if(is.na(mymasks[[1]][1])) {
                    mymasks[[3+1]] <- seq(lowlim, uplim)
                } else {
                    mymasks[[3+1]] <- mymasks[[1]] & mymasks[[2]] & mymasks[[3]]
                }
            }


            ## extrapolate epsilon to zero
            fitresult <- try(bootstrap.nlsfit(x=data[[2+iset]]$epsilons[(epslowlim+1):(data[[2+iset]]$metadata$neps-epsuplim)], y=DGamma,
                                bsamples=bs, fn=fitfn, par.guess=par.guess,
                                mask=mymasks[[iz+1]], verbose=F))
            if(!inherits(fitresult, "try-error")) {
            try(plot(fitresult, main=paste(title, "only stat"), xlab="epsilon", ylab="DG/Dq2"))
            try(plotwitherror(x=0, y=fitresult$t0[1], dy=fitresult$se[1], rep=TRUE, col="red"))

            ## save result
            result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=data[[2+iset]]$metadata$w,
                                        nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=fitresult$t0[1], dDGDq2=fitresult$se[1],
                                        chi=fitresult$chisqr / fitresult$dof, p=fitresult$Qval, includesys=F, fitcorr=FALSE))
            resultdat$tsnk[2*parindex-1]       <- tsnk[index]
            resultdat$Nt[2*parindex-1]         <- Nt[index]
            resultdat$w[2*parindex-1]          <- data[[2+iset]]$metadata$w
            resultdat$nerr[2*parindex-1]       <- nerr[index]
            resultdat$iz[2*parindex-1]         <- iz
            resultdat$icomb[2*parindex-1]      <- icomb
            resultdat$th[2*parindex-1]         <- th[index]
            resultdat$includesys[2*parindex-1] <- FALSE
            resultdat$DGDq2[2*parindex-1]      <- fitresult$t0[1]
            resultdat$dDGDq2[2*parindex-1]     <- fitresult$se[1]
            resultdat$dat[, 2*parindex-1]        <- fitresult$t[, 1]
            
            
            } else {
                result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=data[[2+iset]]$metadata$w, nerr=nerr[index],
                                        iz=iz, icomb=icomb, th=th[index], DGDq2=NA, dDGDq2=NA, chi=NA, p=NA, includesys=NA, fitcorr=NA))
            }


            ## extrapolate with total error
            bootsys <- parametric.bootstrap(boot.R=data[[2+iset]][[namesel[1]]]$bootnumber[[1]], x=rep(0, data[[2+iset]]$metadata$neps-epslowlim-epsuplim), dx=dDGamma_sys, seed=12345678)
            fitresult <- try(bootstrap.nlsfit(x=data[[2+iset]]$epsilons[(epslowlim+1):(data[[2+iset]]$metadata$neps-epsuplim)], y=DGamma, bsamples=bs+bootsys,
                                fn=fitfn, par.guess=par.guess, mask=mymasks[[iz+1]]))
            if(!inherits(fitresult, "try-error")) {

            try(plot(fitresult, main=paste(title, "total error"), xlab="epsilon", ylab="DG/Dq2"))
            try(plotwitherror(x=0, y=fitresult$t0[1], dy=fitresult$se[1], rep=TRUE, col="red"))
            result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=data[[2+iset]]$metadata$w,
                                        nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=fitresult$t0[1], dDGDq2=fitresult$se[1],
                                        chi=fitresult$chisqr / fitresult$dof, p=fitresult$Qval, includesys=T, fitcorr=F))
                                        
            resultdat$tsnk[2*parindex]       <- tsnk[index]
            resultdat$Nt[2*parindex]         <- Nt[index]
            resultdat$w[2*parindex]          <- data[[2+iset]]$metadata$w
            resultdat$nerr[2*parindex]       <- nerr[index]
            resultdat$iz[2*parindex]         <- iz
            resultdat$icomb[2*parindex]      <- icomb
            resultdat$th[2*parindex]         <- th[index]
            resultdat$includesys[2*parindex] <- TRUE
            resultdat$DGDq2[2*parindex]      <- fitresult$t0[1]
            resultdat$dDGDq2[2*parindex]     <- fitresult$se[1]
            resultdat$dat[, 2*parindex]        <- fitresult$t[, 1]
            
            } else {
                result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=data[[2+iset]]$metadata$w,
                                            nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=NA, dDGDq2=NA, chi=NA, p=NA,
                                            includesys=NA, fitcorr=NA))
            }
            parindex <- parindex + 1

        }
    }
}
}
result <- result[-1, ]
result <- result[ ! (result$iz==0 & (result$icomb==3 | result$icomb==4)), ]
print(warnings())
write.table(result, paste0(savename, ".csv"), row.names=F, col.names=T)
saveRDS(resultdat, paste0(savename, ".RDS"))
}


