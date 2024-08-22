#~ library("hadron")
source("/home/gross/Documents/heavymesons/scripts/extrapolate_sigma_templates/read_in_binary.R")

determineDGDq2 <- function(resultpath, filenames, tsnk, Nt, th,
        nerr, amin, savename, fitfn, par.guess, isets=c(0), icomb=c(0),
        epsuplim=0, epslowlim=0) {
    ## detect some errors
    stopifnot(length(filenames) == length(tsnk))
    stopifnot(length(filenames) == length(Nt))
    stopifnot(length(filenames) == length(th))
    stopifnot(length(filenames) == length(nerr))


    result <- data.frame(tsnk=NA, Nt=NA, w=NA, nerr=NA, iz=NA, icomb=NA, th = NA, DGDq2=NA, dDGDq2=NA, chi=NA, p=NA, includesys=NA, fitcorr=NA)

    resultdat <- list(tsnk=c(), Nt=c(), w=c(), nerr=c(), iz=c(), icomb=c(), th=c(), includesys = c(), DGDq2 = c(), dDGDq2 = c(), dat = array(NA, dim=c(1000, length(tsnk)*2*length(isets)*length(icomb)*4)))
    parindex <- 1
for (index in seq(1, length(filenames))){
    data <- read_in_DGDq2(filename=filenames[index], write=FALSE, resultpath=resultpath, NDG=maxz)
for (iset in isets) {
    mymasks <- list(c(NA))
    for (iz in c(0, 1, 2, 3)) {
        for(icomb in icomb) {

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


## bootsample = bootsample - (mean - bootsample) * (errtotmean / bootstat - 1)
increasedistancebootstrap <- function(sample, mean, ratio) (sample - (mean - sample) * (ratio - 1))


determineDGDq2_all <- function(resultpath, filenames, tsnk, Nt, th,
        nerr, amin, savename, fitfn, par.guess, isets=c(0), comblist=c(0), zlist=c(0, 1, 2, 3),
        epsuplim=0, epslowlim=0, errors=c("stat"), volumetable = "",
        doplot=TRUE, NDG = 4, maxz=4, bsamples=1000, neps=-1) {

    ## all members of errors have to be one of "stat", "sys", "vol", "tot" and no duplicates
    ## for calculating to we have to calculate sys and vol as well
    for (i in seq(1, length(errors))){
        stopifnot(errors[i] == "stat" || errors[i] == "sys" || errors[i] == "vol" || errors[i] == "tot" )
    }
    stopifnot(length(unique(errors)) == length(errors))

    if ("vol" %in% errors) stopifnot(file.exists(volumetable))
    stopifnot(!("tot" %in% errors && (!("sys" %in% errors) || !("vol" %in% errors))))

    ## detect some errors
    stopifnot(length(filenames) == length(tsnk))
    stopifnot(length(filenames) == length(Nt))
    stopifnot(length(filenames) == length(th))
    stopifnot(length(filenames) == length(nerr))
    
    ## and some probable mistakes
    if(NDG!=4 && NDG!=5) print("Are you sure about your value of NDG?")
    


    result <- data.frame(tsnk=NA, Nt=NA, w=NA, nerr=NA, iz=NA, icomb=NA, th = NA, DGDq2=NA, dDGDq2=NA, chi=NA, p=NA, errtype=NA)

    resultdat <- list(tsnk=c(), Nt=c(), w=c(), nerr=c(), iz=c(), 
    icomb=c(), th=c(), errtype = c(), DGDq2 = c(), dDGDq2 = c(), 
    dat = array(NA, dim=c(bsamples, length(tsnk)*4*length(isets)*length(comblist)*maxz)))
    
    fitlist <- list()
    namefitlist <- seq(1, length(tsnk)*4*length(isets)*length(comblist)*maxz)

    ## read in table with finite volume results

    if ("vol" %in% errors || "tot" %in% errors) finitevolumeall <- read.table(volumetable, header=TRUE)


    parindex <- 1
for (index in seq(1, length(filenames))){
    data <- try(read_in_DGDq2(filename=filenames[index], write=FALSE, resultpath=resultpath, NDG = NDG))
for (iset in isets) {
    mymasks <- list(c(NA))
    for (iz in zlist) {
        for(icomb in comblist) {
            
            if(!inherits(data, "try-error")) {
            if ("vol" %in% errors || "tot" %in% errors) volume <- finitevolumeall[finitevolumeall$th == th[index] & finitevolumeall$iz == iz, ]
            
            ## set number of epsilons: default epslowlim:epsuplim, if neps is given, epslowlim:(epslowlim+neps)
            stopifnot(data[[2+iset]]$metadata$neps >= epslowlim+neps)
            if(neps!=-1) epsuplim <- data[[2+iset]]$metadata$neps - (epslowlim+neps)

            ## determine which part of data to use, extract needed elements
            namesel <- paste0("id", ((iz)*5 + icomb)*data[[2+iset]]$metadata$neps + (epslowlim:(data[[2+iset]]$metadata$neps-1-epsuplim)), "ieps", epslowlim:(data[[2+iset]]$metadata$neps-1-epsuplim), "icomb", icomb, "iz", iz)
            title=paste("iset", iset, "iz", iz, "icomb", icomb, "tsnk", tsnk[index], "Nt", Nt[index], "th", th[index], "nerr", nerr[index])
            print(title)


            DGamma <- unlist(sapply(X=data[[2+iset]][namesel], FUN=getElement, name="DGDq2mean")[1, ], use.names=F)
            dDGamma_sys <- unname(unlist(sapply(X=data[[2+iset]][namesel], FUN=getElement, name="sys")[1, ], use.names=F))
            dDGamma <- unlist(sapply(X=data[[2+iset]][namesel], FUN=getElement, name="DGDq2sd")[1, ], use.names=F)
            DGammaboot <- sapply(X=data[[2+iset]][namesel], FUN=getElement, name="DGDq2boots")[1, ]
            
            bs <- array(unlist(DGammaboot, use.names=F), dim=c(data[[2+iset]][[namesel[1]]]$bootnumber[[1]], data[[2+iset]]$metadata$neps-epslowlim-epsuplim))
            x <- data[[2+iset]]$epsilons[(epslowlim+1):(data[[2+iset]]$metadata$neps-epsuplim)]

            ## determine min A/A0_ref
            AA0ar <- array(rep(NA, data[[2+iset]]$metadata$neps*data[[2+iset]]$metadata$nnorm), dim=c(data[[2+iset]]$metadata$nnorm, data[[2+iset]]$metadata$neps-epslowlim-epsuplim))
            for (norm in seq(1, data[[2+iset]]$metadata$nnorm)) {
                A0ABCW_ref <- sapply(X=data[[2+iset]][namesel], FUN=getElement, name="A0ABCW_ref")[1, ]
                A0ABCW_ref_ar <- t(array(unlist(A0ABCW_ref, use.names=F), dim=c(6, data[[2+iset]]$metadata$neps-epslowlim-epsuplim)))
                AA0ar[norm, ] <- A0ABCW_ref_ar[, 2] / A0ABCW_ref_ar[, 1]
            }
            AA0 <- apply(X=AA0ar, MARGIN=2, FUN=min)

            ## determine masks for fit
            if(iz != maxz-1) mymasks[[iz+1]] <- (AA0 < amin)
            if(iz == maxz-1) {
                if(is.na(mymasks[[1]][1])) {
                    mymasks[[maxz]] <- seq(epslowlim, data[[2+iset]]$metadata$nnorm - epsuplim)
                } else {
                    mymasks[[maxz]] <- mymasks[[1]] 
                    for (zmask in seq(2, maxz)) {
                        mymasks[[maxz]] <- mymasks[[maxz]] & mymasks[[zmask]]
                    }
                }
            }

            ## error only from statistics
            if ("stat" %in% errors) {
                ## extrapolate epsilon to zero
                
                fitresult <- try(bootstrap.nlsfit(x=data[[2+iset]]$epsilons[(epslowlim+1):(data[[2+iset]]$metadata$neps-epsuplim)], y=DGamma,
                                    bsamples=bs, fn=fitfn, par.guess=par.guess,
                                    mask=mymasks[[iz+1]], verbose=(iz==3)))
                if(!inherits(fitresult, "try-error")) {
                if (doplot) try(plot(fitresult, main=paste(title, "only stat"), xlab="epsilon", ylab="DG/Dq2", xlim=c(0, max(fitresult$x))))
                if (doplot) try(plotwitherror(x=0, y=fitresult$t0[1], dy=fitresult$se[1], rep=TRUE, col="red"))

                ## save result
                result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=data[[2+iset]]$metadata$w,
                                            nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=fitresult$t0[1], dDGDq2=fitresult$se[1],
                                            chi=fitresult$chisqr / fitresult$dof, p=fitresult$Qval, errtype = "stat"))
                resultdat$tsnk[4*parindex-3]       <- tsnk[index]
                resultdat$Nt[4*parindex-3]         <- Nt[index]
                resultdat$w[4*parindex-3]          <- data[[2+iset]]$metadata$w
                resultdat$nerr[4*parindex-3]       <- nerr[index]
                resultdat$iz[4*parindex-3]         <- iz
                resultdat$icomb[4*parindex-3]      <- icomb
                resultdat$th[4*parindex-3]         <- th[index]
                resultdat$errtype[4*parindex-3]    <- "stat"
                resultdat$DGDq2[4*parindex-3]      <- fitresult$t0[1]
                resultdat$dDGDq2[4*parindex-3]     <- fitresult$se[1]
                resultdat$dat[, 4*parindex-3]      <- fitresult$t[, 1]
                fitlist[[4*parindex-3]]            <- fitresult
                namefitlist[4*parindex-3]          <- paste(title, "stat")


                } else {
                    result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=data[[2+iset]]$metadata$w, nerr=nerr[index],
                                            iz=iz, icomb=icomb, th=th[index], DGDq2=NA, dDGDq2=NA, chi=NA, p=NA, errtype=NA))
                    fitlist[[4*parindex-3]] <- NA
                }
            } else fitlist[[4*parindex-3]] <- NA


            ## extrapolate with total error
            if ("sys" %in% errors) {
                ratio <- sqrt(dDGamma^2 + dDGamma_sys^2)/dDGamma
                bootsys <- sweep(x=bs, STAT=-ratio, MARGIN=2, FUN='*') + array(rep(DGamma*(1+ratio), each=bsamples), dim=c(bsamples, length(ratio)))
                        
#~                 print(dDGamma_sys/dDGamma)
                
                fitresult <- try(bootstrap.nlsfit(x=data[[2+iset]]$epsilons[(epslowlim+1):(data[[2+iset]]$metadata$neps-epsuplim)], y=DGamma, bsamples=bootsys,
                                    fn=fitfn, par.guess=par.guess, mask=mymasks[[iz+1]]))
                if(!inherits(fitresult, "try-error")) {

                if (doplot) try(plot(fitresult, main=paste(title, "stat+sys error"), xlab="epsilon", ylab="DG/Dq2", xlim=c(0, max(fitresult$x))))
                if (doplot) try(plotwitherror(x=0, y=fitresult$t0[1], dy=fitresult$se[1], rep=TRUE, col="red"))
                result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=data[[2+iset]]$metadata$w,
                                            nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=fitresult$t0[1], dDGDq2=fitresult$se[1],
                                            chi=fitresult$chisqr / fitresult$dof, p=fitresult$Qval, errtype="sys"))

                resultdat$tsnk[4*parindex-2]       <- tsnk[index]
                resultdat$Nt[4*parindex-2]         <- Nt[index]
                resultdat$w[4*parindex-2]          <- data[[2+iset]]$metadata$w
                resultdat$nerr[4*parindex-2]       <- nerr[index]
                resultdat$iz[4*parindex-2]         <- iz
                resultdat$icomb[4*parindex-2]      <- icomb
                resultdat$th[4*parindex-2]         <- th[index]
                resultdat$errtype[4*parindex-2]    <- "sys"
                resultdat$DGDq2[4*parindex-2]      <- fitresult$t0[1]
                resultdat$dDGDq2[4*parindex-2]     <- fitresult$se[1]
                resultdat$dat[, 4*parindex-2]      <- fitresult$t[, 1]
                fitlist[[4*parindex-2]]            <- fitresult
                namefitlist[4*parindex-2]          <- paste(title, "sys")

                } else {
                    result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=data[[2+iset]]$metadata$w,
                                                nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=NA, dDGDq2=NA, chi=NA, p=NA,
                                                errtype=NA))
                    fitlist[[4*parindex-2]] <- NA
                }
            } else fitlist[[4*parindex-2]] <- NA

            ## extrapolate with finite volume error
            if ("vol" %in% errors) {
                stopifnot(data[[2+iset]]$epsilons == volume$epsilons)
                
                ratio <- sqrt(dDGamma^2 + volume$Delta[(epslowlim+1):(data[[2+iset]]$metadata$neps-epsuplim)]^2)/dDGamma
                bootvol <- sweep(x=bs, STAT=-ratio, MARGIN=2, FUN='*') + array(rep(DGamma*(1+ratio), each=bsamples), dim=c(bsamples, length(ratio)))
                
                fitresult <- try(bootstrap.nlsfit(x=data[[2+iset]]$epsilons[(epslowlim+1):(data[[2+iset]]$metadata$neps-epsuplim)], y=DGamma, bsamples=bootvol,
                                    fn=fitfn, par.guess=par.guess, mask=mymasks[[iz+1]]))
                if(!inherits(fitresult, "try-error")) {

                if (doplot) try(plot(fitresult, main=paste(title, "stat+vol error"), xlab="epsilon", ylab="DG/Dq2", xlim=c(0, max(fitresult$x))))
                if (doplot) try(plotwitherror(x=0, y=fitresult$t0[1], dy=fitresult$se[1], rep=TRUE, col="red"))
                result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=data[[2+iset]]$metadata$w,
                                            nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=fitresult$t0[1], dDGDq2=fitresult$se[1],
                                            chi=fitresult$chisqr / fitresult$dof, p=fitresult$Qval, errtype="vol"))

                resultdat$tsnk[4*parindex-1]       <- tsnk[index]
                resultdat$Nt[4*parindex-1]         <- Nt[index]
                resultdat$w[4*parindex-1]          <- data[[2+iset]]$metadata$w
                resultdat$nerr[4*parindex-1]       <- nerr[index]
                resultdat$iz[4*parindex-1]         <- iz
                resultdat$icomb[4*parindex-1]      <- icomb
                resultdat$th[4*parindex-1]         <- th[index]
                resultdat$errtype[4*parindex-1]    <- "vol"
                resultdat$DGDq2[4*parindex-1]      <- fitresult$t0[1]
                resultdat$dDGDq2[4*parindex-1]     <- fitresult$se[1]
                resultdat$dat[, 4*parindex-1]      <- fitresult$t[, 1]
                fitlist[[4*parindex-1]]            <- fitresult
                namefitlist[4*parindex-1]          <- paste(title, "vol")

                } else {
                    result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=data[[2+iset]]$metadata$w,
                                                nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=NA, dDGDq2=NA, chi=NA, p=NA,
                                                errtype=NA))
                    fitlist[[4*parindex-2]] <- NA
                }
            } else fitlist[[4*parindex-1]] <- NA

            ## extrapolate with finite volume and systematic error
            if ("tot" %in% errors) {
                stopifnot(data[[2+iset]]$epsilons == volume$epsilons)
                
                ratio <- sqrt(dDGamma^2 + dDGamma_sys^2 + volume$Delta[(epslowlim+1):(data[[2+iset]]$metadata$neps-epsuplim)]^2)/dDGamma
                boottot <- sweep(x=bs, STAT=-ratio, MARGIN=2, FUN='*') + array(rep(DGamma*(1+ratio), each=bsamples), dim=c(bsamples, length(ratio)))
                
                fitresult <- try(bootstrap.nlsfit(x=data[[2+iset]]$epsilons[(epslowlim+1):(data[[2+iset]]$metadata$neps-epsuplim)], y=DGamma, bsamples=boottot,
                                    fn=fitfn, par.guess=par.guess, mask=mymasks[[iz+1]]))
                if(!inherits(fitresult, "try-error")) {

                if (doplot) try(plot(fitresult, main=paste(title, "total error"), xlab="epsilon", ylab="DG/Dq2", xlim=c(0, max(fitresult$x))))
                if (doplot) try(plotwitherror(x=0, y=fitresult$t0[1], dy=fitresult$se[1], rep=TRUE, col="red"))
                result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=data[[2+iset]]$metadata$w,
                                            nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=fitresult$t0[1], dDGDq2=fitresult$se[1],
                                            chi=fitresult$chisqr / fitresult$dof, p=fitresult$Qval, errtype="tot"))

                resultdat$tsnk[4*parindex]       <- tsnk[index]
                resultdat$Nt[4*parindex]         <- Nt[index]
                resultdat$w[4*parindex]          <- data[[2+iset]]$metadata$w
                resultdat$nerr[4*parindex]       <- nerr[index]
                resultdat$iz[4*parindex]         <- iz
                resultdat$icomb[4*parindex]      <- icomb
                resultdat$th[4*parindex]         <- th[index]
                resultdat$errtype[4*parindex]    <- "tot"
                resultdat$DGDq2[4*parindex]      <- fitresult$t0[1]
                resultdat$dDGDq2[4*parindex]     <- fitresult$se[1]
                resultdat$dat[, 4*parindex]      <- fitresult$t[, 1]
                fitlist[[4*parindex]]            <- fitresult
                namefitlist[4*parindex]          <- paste(title, "tot")

                } else {
                    result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=data[[2+iset]]$metadata$w,
                                                nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=NA, dDGDq2=NA, chi=NA, p=NA,
                                                errtype=NA))
                    fitlist[[4*parindex]]            <- NA
                    namefitlist[4*parindex]          <- paste(title, "tot")
                }
            } else {
                #we set these elements to NA in the case they are not calculated to ensure all lists are filled completetly
                resultdat$tsnk[4*parindex]       <- NA
                resultdat$Nt[4*parindex]         <- NA
                resultdat$w[4*parindex]          <- NA
                resultdat$nerr[4*parindex]       <- NA
                resultdat$iz[4*parindex]         <- NA
                resultdat$icomb[4*parindex]      <- NA
                resultdat$th[4*parindex]         <- NA
                resultdat$errtype[4*parindex]    <- NA
                resultdat$DGDq2[4*parindex]      <- NA
                resultdat$dDGDq2[4*parindex]     <- NA
                resultdat$dat[, 4*parindex]      <- rep(NA, bsamples)
                fitlist[[4*parindex]]            <- NA
                namefitlist[4*parindex]          <- NA
            }



            parindex <- parindex + 1

    } else {
                #we set these elements to NA in the case they are not calculated to ensure all lists are filled completetly
                resultdat$tsnk[4*parindex]       <- NA
                resultdat$Nt[4*parindex]         <- NA
                resultdat$w[4*parindex]          <- NA
                resultdat$nerr[4*parindex]       <- NA
                resultdat$iz[4*parindex]         <- NA
                resultdat$icomb[4*parindex]      <- NA
                resultdat$th[4*parindex]         <- NA
                resultdat$errtype[4*parindex]    <- NA
                resultdat$DGDq2[4*parindex]      <- NA
                resultdat$dDGDq2[4*parindex]     <- NA
                resultdat$dat[, 4*parindex]      <- rep(NA, bsamples)
                fitlist[[4*parindex]]            <- NA
                namefitlist[4*parindex]          <- NA
            parindex <- parindex + 1
            }
        }
    }
}
}
result <- result[-1, ]
result <- result[ ! (result$iz==0 & (result$icomb==3 | result$icomb==4)), ]
print(warnings())
write.table(result, paste0(savename, ".csv"), row.names=F, col.names=T)
saveRDS(resultdat, paste0(savename, ".RDS"))

names(fitlist) <- namefitlist
return(invisible(fitlist))
}

