
source("/home/gross/Documents/heavymesons/scripts_christiane/extrapolate_sigma_templates/read_in_binary.R")
source("/home/gross/Documents/heavymesons/scripts_christiane/extrapolate_sigma_templates/aichelpers.R")

## bootsample = bootsample - (mean - bootsample) * (errtotmean / bootstat - 1)
increasedistancebootstrap <- function(sample, mean, ratio) (sample - (mean - sample) * (ratio - 1))


## TODO: different finite volume effects possible
determineDGDq2_combine_aic <- function(resultpathlist, filenames, tsnk, Nt, th,
        nerr, amin, savename, fitfnlist, par.guesslist, isets=c(0), comblist=c(0), zlist=c(0, 1, 2, 3),
        epsuplim=0, epslowlim=0, errors=c("stat"), volumetable = "",
        doplot=TRUE, NDG = 4, maxz=4, bsamples=1000, numberfits=2) {

    ## all members of errors have to be one of "stat", "sys", "vol", "tot" and no duplicates
    ## for calculating tot we have to calculate sys and vol as well
    for (i in seq(1, length(errors))){
        stopifnot(errors[i] == "stat" || errors[i] == "sys" || errors[i] == "vol" || errors[i] == "tot" )
    }

    if ("vol" %in% errors) stopifnot(file.exists(volumetable))
    stopifnot(!("tot" %in% errors && (!("sys" %in% errors) || !("vol" %in% errors))))

    ## detect some errors
    stopifnot(length(unique(errors)) == length(errors))
    stopifnot(numberfits==length(fitfnlist))
    stopifnot(numberfits==length(par.guesslist))
    stopifnot(numberfits==length(resultpathlist))
    stopifnot(length(filenames)%%numberfits == 0)
    stopifnot(length(filenames)/numberfits == length(tsnk))
    stopifnot(length(filenames)/numberfits == length(Nt))
    stopifnot(length(filenames)/numberfits == length(th))
    stopifnot(length(filenames)/numberfits == length(nerr))
    
    ## and some probable mistakes
    if(NDG!=4 && NDG!=5) print("Are you sure about your value of NDG?")


    result <- data.frame(tsnk=NA, Nt=NA, w=NA, nerr=NA, iz=NA, icomb=NA, th = NA, DGDq2=NA, dDGDq2=NA, errtype=NA)

    resultdat <- list(tsnk=c(), Nt=c(), w=c(), nerr=c(), iz=c(), 
    icomb=c(), th=c(), errtype = c(), DGDq2 = c(), dDGDq2 = c(), 
    dat = array(NA, dim=c(bsamples, length(tsnk)*4*length(isets)*length(comblist)*maxz)))
    
    fitlist <- list()
    namefitlist <- seq(1, length(tsnk)*4*length(isets)*length(comblist)*maxz)

    ## read in table with finite volume results

    if ("vol" %in% errors || "tot" %in% errors) finitevolumeall <- read.table(volumetable, header=TRUE)


    parindex <- 1
    datalist <- list()
for (index in seq(1, length(filenames)/numberfits)){
    for (datanumber in seq(1, numberfits)) {
        print(filenames[index+(datanumber-1)*length(tsnk)])
        datalist[[datanumber]] <- try(read_in_DGDq2(filename=filenames[index+(datanumber-1)*length(tsnk)], write=FALSE, resultpath=resultpathlist[[datanumber]], NDG = NDG))
    }
for (iset in isets) {
    mymasks <- list(c(NA))
    for (iz in zlist) {
        for(icomb in comblist) {
            do_continue <- TRUE
            for (datanumber in seq(1, numberfits)) {
                do_continue <- do_continue && !inherits(data, "try-error")
            }
            if(do_continue) {
            if ("vol" %in% errors || "tot" %in% errors) volume <- finitevolumeall[finitevolumeall$th == th[index] & finitevolumeall$iz == iz, ]

            DGammalist <- list()
            dDGammasyslist <- list()
            bslist <- list()
            xlist <- list()
            

            for (datanumber in seq(1, numberfits)) {
            ## determine which part of data to use, extract needed elements
            namesel <- paste0("id", ((iz)*5 + icomb)*datalist[[datanumber]][[2+iset]]$metadata$neps + (epslowlim:(datalist[[datanumber]][[2+iset]]$metadata$neps-1-epsuplim)), "ieps", epslowlim:(datalist[[datanumber]][[2+iset]]$metadata$neps-1-epsuplim), "icomb", icomb, "iz", iz)
            title=paste("iset", iset, "iz", iz, "icomb", icomb, "tsnk", tsnk[index], "Nt", Nt[index], "th", th[index], "nerr", nerr[index])
            print(title)

## TODO: make lists for everything, apply AIC
            DGammalist[[datanumber]] <- unlist(sapply(X=datalist[[datanumber]][[2+iset]][namesel], FUN=getElement, name="DGDq2mean")[1, ], use.names=F)
            dDGammasyslist[[datanumber]] <- unname(unlist(sapply(X=datalist[[datanumber]][[2+iset]][namesel], FUN=getElement, name="sys")[1, ], use.names=F))
            dDGamma <- sqrt(unlist(sapply(X=datalist[[datanumber]][[2+iset]][namesel], FUN=getElement, name="DGDq2sd")[1, ], use.names=F))
            DGammaboot <- sapply(X=datalist[[datanumber]][[2+iset]][namesel], FUN=getElement, name="DGDq2boots")[1, ]
            
            bslist[[datanumber]] <- array(unlist(DGammaboot, use.names=F), dim=c(datalist[[datanumber]][[2+iset]][[namesel[1]]]$bootnumber[[1]], datalist[[datanumber]][[2+iset]]$metadata$neps-epslowlim-epsuplim))
            xlist[[datanumber]] <- datalist[[datanumber]][[2+iset]]$epsilons[(epslowlim+1):(datalist[[datanumber]][[2+iset]]$metadata$neps-epsuplim)]
            
            ## determine min A/A0_ref
            AA0ar <- array(NA, dim=c(datalist[[datanumber]][[2+iset]]$metadata$nnorm, datalist[[datanumber]][[2+iset]]$metadata$neps-epslowlim-epsuplim))
            for (norm in seq(1, datalist[[datanumber]][[2+iset]]$metadata$nnorm)) {
                A0ABCW_ref <- sapply(X=datalist[[datanumber]][[2+iset]][namesel], FUN=getElement, name="A0ABCW_ref")[1, ]
                A0ABCW_ref_ar <- t(array(unlist(A0ABCW_ref, use.names=F), dim=c(6, datalist[[datanumber]][[2+iset]]$metadata$neps-epslowlim-epsuplim)))
                AA0ar[norm, ] <- A0ABCW_ref_ar[, 2] / A0ABCW_ref_ar[, 1]
            }
            AA0 <- apply(X=AA0ar, MARGIN=2, FUN=min)

            ## determine masks for fit
            ## how to do this with lists?
            ## fix this, indices are wrong!
#~             print((datanumber-1)*maxz+iz+1)
            if(iz != maxz-1) mymasks[[(datanumber-1)*maxz+iz+1]] <- (AA0 < amin)
            if (iz==maxz-1) {
                if(is.na(mymasks[[(datanumber-1)*maxz+0+1]][1])) {
                    mymasks[[(datanumber-1)*maxz+iz+1]] <- seq(epslowlim, data[[2+iset]]$metadata$nnorm - epsuplim)
                } else {
                    mymasks[[(datanumber-1)*maxz+iz+1]] <- mymasks[[(datanumber-1)*maxz+0+1]]
                    for (zmask in seq(1, maxz-2)) {
                        mymasks[[(datanumber-1)*maxz+iz+1]] <- mymasks[[(datanumber-1)*maxz+iz+1]] & mymasks[[(datanumber-1)*maxz+zmask+1]]
                    }
                }
            }
        }

            ## error only from statistics
            if ("stat" %in% errors) {
                ## extrapolate epsilon to zero
                fitresult <- try(combinefunctionsaic(xlist=xlist, y=DGammalist,
                                    bsamples=bslist, functions=fitfnlist, pars=par.guesslist,
                                    mask=mymasks[(0:(numberfits-1))*maxz+iz+1], verbose=T))
                
                if(!inherits(fitresult, "try-error")) {
                if (doplot) try(plotlistfits(fitresult$fitresults, main=paste(title, "only stat"), xlab="epsilon", ylab="DG/Dq2", cols=c("blue", "green")))
                if (doplot) try(plotwitherror(x=0, y=fitresult$cdfsummary$t0[1], dy=fitresult$cdfsummary$errstat[1], rep=TRUE, col="red", cex=2, lwd=2))

                ## save result
                result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=data[[2+iset]]$metadata$w,
                                            nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=fitresult$cdfsummary$t0[1], dDGDq2=fitresult$cdfsummary$errstat[1],
                                            errtype = "stat"))
                resultdat$tsnk[4*parindex-3]       <- tsnk[index]
                resultdat$Nt[4*parindex-3]         <- Nt[index]
                resultdat$w[4*parindex-3]          <- datalist[[1]][[2+iset]]$metadata$w
                resultdat$nerr[4*parindex-3]       <- nerr[index]
                resultdat$iz[4*parindex-3]         <- iz
                resultdat$icomb[4*parindex-3]      <- icomb
                resultdat$th[4*parindex-3]         <- th[index]
                resultdat$errtype[4*parindex-3]    <- "stat"
                resultdat$DGDq2[4*parindex-3]      <- fitresult$cdfsummary$t0[1]
                resultdat$dDGDq2[4*parindex-3]     <- fitresult$cdfsummary$errstat[1]
                resultdat$dat[, 4*parindex-3]      <- fitresult$cdfbootsamples$bootmean[, 1]
                fitlist[[4*parindex-3]]            <- fitresult
                namefitlist[4*parindex-3]          <- paste(title, "stat")


#~                 } else {
#~                     result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=data[[2+iset]]$metadata$w, nerr=nerr[index],
#~                                             iz=iz, icomb=icomb, th=th[index], DGDq2=NA, dDGDq2=NA, chi=NA, p=NA, errtype=NA))
#~                     fitlist[[4*parindex-3]] <- NA
                }
            } else fitlist[[4*parindex-3]] <- NA


#~             ## extrapolate with total error
#~             if ("sys" %in% errors) {
#                bootsys <- parametric.bootstrap(boot.R=data[[2+iset]][[namesel[1]]]$bootnumber[[1]], 
#                        x=rep(0, data[[2+iset]]$metadata$neps-epslowlim-epsuplim), 
#                        dx=dDGamma_sys, seed=12345678)
#~ ## bootsample = bootsample - (mean - bootsample) * (errtotmean / bootstat - 1)
#~                 bootsys <- t(apply(X=bs, MARGIN=1, FUN=increasedistancebootstrap, mean=DGamma, 
#~                         ratio=sqrt(dDGamma^2 + dDGamma_sys^2)/dDGamma))
                        
#                print(dDGamma_sys/dDGamma)
                
#~                 fitresult <- try(bootstrap.nlsfit(x=data[[2+iset]]$epsilons[(epslowlim+1):(data[[2+iset]]$metadata$neps-epsuplim)], y=DGamma, bsamples=bootsys,
#~                                     fn=fitfn, par.guess=par.guess, mask=mymasks[[iz+1]]))
#~                 if(!inherits(fitresult, "try-error")) {

#~                 if (doplot) try(plot(fitresult, main=paste(title, "stat+sys error"), xlab="epsilon", ylab="DG/Dq2"))
#~                 if (doplot) try(plotwitherror(x=0, y=fitresult$t0[1], dy=fitresult$se[1], rep=TRUE, col="red"))
#~                 result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=data[[2+iset]]$metadata$w,
#~                                             nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=fitresult$t0[1], dDGDq2=fitresult$se[1],
#~                                             chi=fitresult$chisqr / fitresult$dof, p=fitresult$Qval, errtype="sys"))

#~                 resultdat$tsnk[4*parindex-2]       <- tsnk[index]
#~                 resultdat$Nt[4*parindex-2]         <- Nt[index]
#~                 resultdat$w[4*parindex-2]          <- data[[2+iset]]$metadata$w
#~                 resultdat$nerr[4*parindex-2]       <- nerr[index]
#~                 resultdat$iz[4*parindex-2]         <- iz
#~                 resultdat$icomb[4*parindex-2]      <- icomb
#~                 resultdat$th[4*parindex-2]         <- th[index]
#~                 resultdat$errtype[4*parindex-2]    <- "sys"
#~                 resultdat$DGDq2[4*parindex-2]      <- fitresult$t0[1]
#~                 resultdat$dDGDq2[4*parindex-2]     <- fitresult$se[1]
#~                 resultdat$dat[, 4*parindex-2]      <- fitresult$t[, 1]
#~                 fitlist[[4*parindex-2]]            <- fitresult
#~                 namefitlist[4*parindex-2]          <- paste(title, "sys")

#~                 } else {
#~                     result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=data[[2+iset]]$metadata$w,
#~                                                 nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=NA, dDGDq2=NA, chi=NA, p=NA,
#~                                                 errtype=NA))
#~                     fitlist[[4*parindex-2]] <- NA
#~                 }
#~             } else fitlist[[4*parindex-2]] <- NA

#~             ## extrapolate with finite volume error
#~             if ("vol" %in% errors) {
#~                 stopifnot(data[[2+iset]]$epsilons == volume$epsilons)
                
#                bootvol <- parametric.bootstrap(boot.R=data[[2+iset]][[namesel[1]]]$bootnumber[[1]], 
#                        x=rep(0, data[[2+iset]]$metadata$neps-epslowlim-epsuplim), 
#                        dx=volume$Delta[epslowlim:(data[[2+iset]]$metadata$neps-epsuplim)], seed=123456789)

#~                 bootvol <- t(apply(X=bs, MARGIN=1, FUN=increasedistancebootstrap, mean=DGamma, 
#~                         ratio=sqrt(dDGamma^2 + volume$Delta[(epslowlim+1):(data[[2+iset]]$metadata$neps-epsuplim)]^2)/dDGamma))
                
#~                 fitresult <- try(bootstrap.nlsfit(x=data[[2+iset]]$epsilons[(epslowlim+1):(data[[2+iset]]$metadata$neps-epsuplim)], y=DGamma, bsamples=bootvol,
#~                                     fn=fitfn, par.guess=par.guess, mask=mymasks[[iz+1]]))
#~                 if(!inherits(fitresult, "try-error")) {

#~                 if (doplot) try(plot(fitresult, main=paste(title, "stat+vol error"), xlab="epsilon", ylab="DG/Dq2"))
#~                 if (doplot) try(plotwitherror(x=0, y=fitresult$t0[1], dy=fitresult$se[1], rep=TRUE, col="red"))
#~                 result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=data[[2+iset]]$metadata$w,
#~                                             nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=fitresult$t0[1], dDGDq2=fitresult$se[1],
#~                                             chi=fitresult$chisqr / fitresult$dof, p=fitresult$Qval, errtype="vol"))

#~                 resultdat$tsnk[4*parindex-1]       <- tsnk[index]
#~                 resultdat$Nt[4*parindex-1]         <- Nt[index]
#~                 resultdat$w[4*parindex-1]          <- data[[2+iset]]$metadata$w
#~                 resultdat$nerr[4*parindex-1]       <- nerr[index]
#~                 resultdat$iz[4*parindex-1]         <- iz
#~                 resultdat$icomb[4*parindex-1]      <- icomb
#~                 resultdat$th[4*parindex-1]         <- th[index]
#~                 resultdat$errtype[4*parindex-1]    <- "vol"
#~                 resultdat$DGDq2[4*parindex-1]      <- fitresult$t0[1]
#~                 resultdat$dDGDq2[4*parindex-1]     <- fitresult$se[1]
#~                 resultdat$dat[, 4*parindex-1]      <- fitresult$t[, 1]
#~                 fitlist[[4*parindex-1]]            <- fitresult
#~                 namefitlist[4*parindex-1]          <- paste(title, "vol")

#~                 } else {
#~                     result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=data[[2+iset]]$metadata$w,
#~                                                 nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=NA, dDGDq2=NA, chi=NA, p=NA,
#~                                                 errtype=NA))
#~                     fitlist[[4*parindex-2]] <- NA
#~                 }
#~             } else fitlist[[4*parindex-1]] <- NA

#~             ## extrapolate with finite volume and systematic error
#~             if ("tot" %in% errors) {
#~                 stopifnot(data[[2+iset]]$epsilons == volume$epsilons)
                
#~                 boottot <- t(apply(X=bs, MARGIN=1, FUN=increasedistancebootstrap, mean=DGamma, 
#~                             ratio=sqrt(dDGamma^2 + dDGamma_sys^2 + volume$Delta[(epslowlim+1):(data[[2+iset]]$metadata$neps-epsuplim)]^2)/dDGamma))
                
#~                 fitresult <- try(bootstrap.nlsfit(x=data[[2+iset]]$epsilons[(epslowlim+1):(data[[2+iset]]$metadata$neps-epsuplim)], y=DGamma, bsamples=boottot,
#~                                     fn=fitfn, par.guess=par.guess, mask=mymasks[[iz+1]]))
#~                 if(!inherits(fitresult, "try-error")) {

#~                 if (doplot) try(plot(fitresult, main=paste(title, "total error"), xlab="epsilon", ylab="DG/Dq2"))
#~                 if (doplot) try(plotwitherror(x=0, y=fitresult$t0[1], dy=fitresult$se[1], rep=TRUE, col="red"))
#~                 result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=data[[2+iset]]$metadata$w,
#~                                             nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=fitresult$t0[1], dDGDq2=fitresult$se[1],
#~                                             chi=fitresult$chisqr / fitresult$dof, p=fitresult$Qval, errtype="tot"))

#~                 resultdat$tsnk[4*parindex]       <- tsnk[index]
#~                 resultdat$Nt[4*parindex]         <- Nt[index]
#~                 resultdat$w[4*parindex]          <- data[[2+iset]]$metadata$w
#~                 resultdat$nerr[4*parindex]       <- nerr[index]
#~                 resultdat$iz[4*parindex]         <- iz
#~                 resultdat$icomb[4*parindex]      <- icomb
#~                 resultdat$th[4*parindex]         <- th[index]
#~                 resultdat$errtype[4*parindex]    <- "tot"
#~                 resultdat$DGDq2[4*parindex]      <- fitresult$t0[1]
#~                 resultdat$dDGDq2[4*parindex]     <- fitresult$se[1]
#~                 resultdat$dat[, 4*parindex]      <- fitresult$t[, 1]
#~                 fitlist[[4*parindex]]            <- fitresult
#~                 namefitlist[4*parindex]          <- paste(title, "tot")

#~                 } else {
#~                     result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=data[[2+iset]]$metadata$w,
#~                                                 nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=NA, dDGDq2=NA, chi=NA, p=NA,
#~                                                 errtype=NA))
#~                     fitlist[[4*parindex]]            <- NA
#~                     namefitlist[4*parindex]          <- paste(title, "tot")
#~                 }
#~             } else {
#~                 #we set these elements to NA in the case they are not calculated to ensure all lists are filled completetly
#~                 resultdat$tsnk[4*parindex]       <- NA
#~                 resultdat$Nt[4*parindex]         <- NA
#~                 resultdat$w[4*parindex]          <- NA
#~                 resultdat$nerr[4*parindex]       <- NA
#~                 resultdat$iz[4*parindex]         <- NA
#~                 resultdat$icomb[4*parindex]      <- NA
#~                 resultdat$th[4*parindex]         <- NA
#~                 resultdat$errtype[4*parindex]    <- NA
#~                 resultdat$DGDq2[4*parindex]      <- NA
#~                 resultdat$dDGDq2[4*parindex]     <- NA
#~                 resultdat$dat[, 4*parindex]      <- rep(NA, bsamples)
#~                 fitlist[[4*parindex]]            <- NA
#~                 namefitlist[4*parindex]          <- NA
#~             }



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

