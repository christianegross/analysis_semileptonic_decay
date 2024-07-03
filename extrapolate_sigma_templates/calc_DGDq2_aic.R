
source("/home/gross/Documents/heavymesons/scripts_christiane/extrapolate_sigma_templates/read_in_binary.R")
source("/home/gross/Documents/heavymesons/scripts_christiane/extrapolate_sigma_templates/aichelpers.R")

## bootsample = bootsample - (mean - bootsample) * (errtotmean / bootstat - 1)
increasedistancebootstrap <- function(sample, mean, ratio) (sample - (mean - sample) * (ratio - 1))


## TODO: different finite volume effects possible
determineDGDq2_combine_aic <- function(resultpathlist, filenames, tsnk, Nt, th,
        nerr, amin, savename, fitfnlist, par.guesslist, isets=c(0), comblist=c(0), zlist=c(0, 1, 2, 3),
        epsuplim=0, epslowlim=0, errors=c("stat"), volumetablelist = c(""),
        doplot=TRUE, NDG = 4, maxz=4, bsamples=1000, numberfits=2, 
        cols=c("blue", "green", "#D2691E", "#556B2F", "cyan", "pink", "red", rep("black", 100)),
        legendargs=list()) {

    ## all members of errors have to be one of "stat", "sys", "vol", "tot" and no duplicates
    ## for calculating tot we have to calculate sys and vol as well
    for (i in seq(1, length(errors))){
        stopifnot(errors[i] == "stat" || errors[i] == "sys" || errors[i] == "vol" || errors[i] == "tot" )
    }

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
    
    finitevolumeall <- list()
    if ("vol" %in% errors || "tot" %in% errors) {
        stopifnot(length(volumetablelist)==numberfits)
        for (i in seq(1, numberfits)) {
            stopifnot(file.exists(volumetablelist[1]))
    ## read in table with finite volume results
            finitevolumeall[[i]] <- read.table(volumetablelist[i], header=TRUE)
        }
    }


    result <- data.frame(tsnk=NA, Nt=NA, w=NA, nerr=NA, iz=NA, icomb=NA, th = NA, DGDq2=NA, dDGDq2=NA, errtype=NA, fitno=NA)

    resultdat <- list(tsnk=c(), Nt=c(), w=c(), nerr=c(), iz=c(), 
    icomb=c(), th=c(), errtype = c(), DGDq2 = c(), dDGDq2 = c(), 
    dat = array(NA, dim=c(bsamples, length(tsnk)*4*length(isets)*length(comblist)*maxz)))
    
    fitlist <- list()
    namefitlist <- seq(1, length(tsnk)*4*length(isets)*length(comblist)*maxz)
    
    legendargs[["col"]]=cols[1:numberfits]
    legendargs[["pch"]]=rep(1, numberfits)



    parindex <- 1
    datalist <- list()
for (index in seq(1, length(filenames)/numberfits)){
    for (datanumber in seq(1, numberfits)) {
#~         print(filenames[index+(datanumber-1)*length(tsnk)])
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

            DGammalist <- list()
            dDGammalist <- list()
            dDGammasyslist <- list()
            bslist <- list()
            xlist <- list()
            bootsyslist <- list()
            bootvollist <- list()
            boottotlist <- list()
            

            for (datanumber in seq(1, numberfits)) {
            ## determine which part of data to use, extract needed elements
            namesel <- paste0("id", ((iz)*5 + icomb)*datalist[[datanumber]][[2+iset]]$metadata$neps + (epslowlim:(datalist[[datanumber]][[2+iset]]$metadata$neps-1-epsuplim)), "ieps", epslowlim:(datalist[[datanumber]][[2+iset]]$metadata$neps-1-epsuplim), "icomb", icomb, "iz", iz)
            title=paste("iset", iset, "iz", iz, "icomb", icomb, "tsnk", tsnk[index], "Nt", Nt[index], "th", th[index], "nerr", nerr[index])
            if(datanumber==1) print(title)

## TODO: make lists for everything, apply AIC
            DGammalist[[datanumber]] <- unlist(sapply(X=datalist[[datanumber]][[2+iset]][namesel], FUN=getElement, name="DGDq2mean")[1, ], use.names=F)
            dDGammasyslist[[datanumber]] <- unname(unlist(sapply(X=datalist[[datanumber]][[2+iset]][namesel], FUN=getElement, name="sys")[1, ], use.names=F))
            dDGammalist[[datanumber]] <- unlist(sapply(X=datalist[[datanumber]][[2+iset]][namesel], FUN=getElement, name="DGDq2sd")[1, ], use.names=F)
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
                                    mask=mymasks[(0:(numberfits-1))*maxz+iz+1], verbose=F))
                
                if(!inherits(fitresult, "try-error")) {
                if (doplot) try(plotlistfits(fitresult$fitresults, main=paste(title, "only stat"), xlab="epsilon", ylab="DG/Dq2", cols=cols, drawlegend=T, legendargs=legendargs))
                if (doplot) try(plotwitherror(x=0, y=fitresult$cdfsummary$t0[1], dy=fitresult$cdfsummary$errstat[1], rep=TRUE, col="red", cex=1, lwd=3))

                ## save result
                result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=datalist[[datanumber]][[2+iset]]$metadata$w,
                                            nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=fitresult$cdfsummary$t0[1], dDGDq2=fitresult$cdfsummary$errstat[1],
                                            errtype = "stat", fitno=-1))
                result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=datalist[[datanumber]][[2+iset]]$metadata$w,
                                            nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=fitresult$savefits$mean, dDGDq2=fitresult$savefits$sd,
                                            errtype = "stat", fitno=fitresult$savefits$index))
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
                resultdat$dat[, 4*parindex-3]      <- fitresult$cdfbootsamples$bootmean
                fitlist[[4*parindex-3]]            <- fitresult
                namefitlist[4*parindex-3]          <- paste(title, "stat")


                } else {
                    result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=datalist[[datanumber]][[2+iset]]$metadata$w, nerr=nerr[index],
                                            iz=iz, icomb=icomb, th=th[index], DGDq2=NA, dDGDq2=NA, chi=NA, p=NA, errtype=NA))
                    fitlist[[4*parindex-3]] <- NA
                }
            } else fitlist[[4*parindex-3]] <- NA


#~             ## extrapolate with systematical error from HLT
            if ("sys" %in% errors) {
            for (datanumber in seq(1, numberfits)) {
                ratio <- sqrt(dDGammalist[[datanumber]]^2 + dDGammasyslist[[datanumber]]^2)/dDGammalist[[datanumber]]
                bootsyslist[[datanumber]] <- sweep(x=bs, STAT=-ratio, MARGIN=2, FUN='*') + array(rep(DGammalist[[datanumber]]*(1+ratio), each=bsamples), dim=c(bsamples, length(ratio)))
                    }
                        
#~ #                print(dDGamma_sys/dDGamma)
                
                fitresult <- try(combinefunctionsaic(xlist=xlist, y=DGammalist,
                                    bsamples=bootsyslist, functions=fitfnlist, pars=par.guesslist,
                                    mask=mymasks[(0:(numberfits-1))*maxz+iz+1], verbose=F))
                if(!inherits(fitresult, "try-error")) {

                if (doplot) try(plotlistfits(fitresult$fitresults, main=paste(title, "stat + sys"), xlab="epsilon", ylab="DG/Dq2", cols=cols, drawlegend=T, legendargs=legendargs))
                if (doplot) try(plotwitherror(x=0, y=fitresult$cdfsummary$t0[1], dy=fitresult$cdfsummary$errstat[1], rep=TRUE, col="red", cex=1, lwd=3))
                result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=datalist[[datanumber]][[2+iset]]$metadata$w,
                                            nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=fitresult$cdfsummary$t0[1], dDGDq2=fitresult$cdfsummary$errstat[1],
                                            errtype = "sys", fitno=-1))
                result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=datalist[[datanumber]][[2+iset]]$metadata$w,
                                            nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=fitresult$savefits$mean, dDGDq2=fitresult$savefits$sd,
                                            errtype = "sys", fitno=fitresult$savefits$index))
                resultdat$tsnk[4*parindex-2]       <- tsnk[index]
                resultdat$Nt[4*parindex-2]         <- Nt[index]
                resultdat$w[4*parindex-2]          <- datalist[[datanumber]][[2+iset]]$metadata$w
                resultdat$nerr[4*parindex-2]       <- nerr[index]
                resultdat$iz[4*parindex-2]         <- iz
                resultdat$icomb[4*parindex-2]      <- icomb
                resultdat$th[4*parindex-2]         <- th[index]
                resultdat$errtype[4*parindex-2]    <- "sys"
                resultdat$DGDq2[4*parindex-2]      <- fitresult$cdfsummary$t0[1]
                resultdat$dDGDq2[4*parindex-2]     <- fitresult$cdfsummary$errstat[1]
                resultdat$dat[, 4*parindex-2]      <- fitresult$cdfbootsamples$bootmean
                fitlist[[4*parindex-2]]            <- fitresult
                namefitlist[4*parindex-2]          <- paste(title, "sys")

                } else {
                    result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=datalist[[datanumber]][[2+iset]]$metadata$w,
                                                nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=NA, dDGDq2=NA, chi=NA, p=NA,
                                                errtype=NA))
                    fitlist[[4*parindex-2]] <- NA
                }
            } else fitlist[[4*parindex-2]] <- NA

            ## extrapolate with finite volume error
            if ("vol" %in% errors) {
                for (datanumber in seq(1, numberfits)) {
                    volume <- finitevolumeall[[datanumber]][finitevolumeall[[datanumber]]$th == th[index] & finitevolumeall[[datanumber]]$iz == iz, ]
                    
                    stopifnot(datalist[[datanumber]][[2+iset]]$epsilons == volume$epsilons)
                    
    
                ratio <- sqrt(dDGammalist[[datanumber]]^2 + volume$Delta[(epslowlim+1):(datalist[[datanumber]][[2+iset]]$metadata$neps-epsuplim)]^2)/dDGammalist[[datanumber]]
                bootvollist[[datanumber]] <- sweep(x=bs, STAT=-ratio, MARGIN=2, FUN='*') + array(rep(DGammalist[[datanumber]]*(1+ratio), each=bsamples), dim=c(bsamples, length(ratio)))
                }
                    
                fitresult <- try(combinefunctionsaic(xlist=xlist, y=DGammalist,
                                    bsamples=bootvollist, functions=fitfnlist, pars=par.guesslist,
                                    mask=mymasks[(0:(numberfits-1))*maxz+iz+1], verbose=F))
                if(!inherits(fitresult, "try-error")) {

                if (doplot) try(plotlistfits(fitresult$fitresults, main=paste(title, "stat + vol"), xlab="epsilon", ylab="DG/Dq2", cols=cols, drawlegend=T, legendargs=legendargs))
                if (doplot) try(plotwitherror(x=0, y=fitresult$cdfsummary$t0[1], dy=fitresult$cdfsummary$errstat[1], rep=TRUE, col="red", cex=1, lwd=3))
                result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=datalist[[datanumber]][[2+iset]]$metadata$w,
                                            nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=fitresult$cdfsummary$t0[1], dDGDq2=fitresult$cdfsummary$errstat[1],
                                            errtype = "vol", fitno=-1))
                result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=datalist[[datanumber]][[2+iset]]$metadata$w,
                                            nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=fitresult$savefits$mean, dDGDq2=fitresult$savefits$sd,
                                            errtype = "vol", fitno=fitresult$savefits$index))

                resultdat$tsnk[4*parindex-1]       <- tsnk[index]
                resultdat$Nt[4*parindex-1]         <- Nt[index]
                resultdat$w[4*parindex-1]          <- datalist[[datanumber]][[2+iset]]$metadata$w
                resultdat$nerr[4*parindex-1]       <- nerr[index]
                resultdat$iz[4*parindex-1]         <- iz
                resultdat$icomb[4*parindex-1]      <- icomb
                resultdat$th[4*parindex-1]         <- th[index]
                resultdat$errtype[4*parindex-1]    <- "vol"
                resultdat$DGDq2[4*parindex-1]      <- fitresult$cdfsummary$t0[1]
                resultdat$dDGDq2[4*parindex-1]     <- fitresult$cdfsummary$errstat[1]
                resultdat$dat[, 4*parindex-1]      <- fitresult$cdfbootsamples$bootmean
                fitlist[[4*parindex-1]]            <- fitresult
                namefitlist[4*parindex-1]          <- paste(title, "vol")

                } else {
                    result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=datalist[[datanumber]][[2+iset]]$metadata$w,
                                                nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=NA, dDGDq2=NA, chi=NA, p=NA,
                                                errtype=NA))
                    fitlist[[4*parindex-2]] <- NA
                }
            } else fitlist[[4*parindex-1]] <- NA

            ## extrapolate with finite volume and systematic error
            if ("tot" %in% errors) {
                for (datanumber in seq(1, numberfits)) {
                    volume <- finitevolumeall[[datanumber]][finitevolumeall[[datanumber]]$th == th[index] & finitevolumeall[[datanumber]]$iz == iz, ]
                    stopifnot(datalist[[datanumber]][[2+iset]]$epsilons == volume$epsilons)
                    
                    ratio <- sqrt(dDGammalist[[datanumber]]^2 + dDGammasyslist[[datanumber]]^2 + volume$Delta[(epslowlim+1):(datalist[[datanumber]][[2+iset]]$metadata$neps-epsuplim)]^2)/dDGammalist[[datanumber]]
                    boottotlist[[datanumber]] <- sweep(x=bs, STAT=-ratio, MARGIN=2, FUN='*') + array(rep(DGammalist[[datanumber]]*(1+ratio), each=bsamples), dim=c(bsamples, length(ratio)))
                }
                
                
                fitresult <- try(combinefunctionsaic(xlist=xlist, y=DGammalist,
                                    bsamples=boottotlist, functions=fitfnlist, pars=par.guesslist,
                                    mask=mymasks[(0:(numberfits-1))*maxz+iz+1], verbose=F))
                if(!inherits(fitresult, "try-error")) {

                if (doplot) try(plotlistfits(fitresult$fitresults, main=paste(title, "stat + vol + sys"), xlab="epsilon", ylab="DG/Dq2", cols=cols, drawlegend=T, legendargs=legendargs))
                if (doplot) try(plotwitherror(x=0, y=fitresult$cdfsummary$t0[1], dy=fitresult$cdfsummary$errstat[1], rep=TRUE, col="red", cex=1, lwd=3))
                result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=datalist[[datanumber]][[2+iset]]$metadata$w,
                                            nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=fitresult$cdfsummary$t0[1], dDGDq2=fitresult$cdfsummary$errstat[1],
                                            errtype = "tot", fitno=-1))
                result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=datalist[[datanumber]][[2+iset]]$metadata$w,
                                            nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], DGDq2=fitresult$savefits$mean, dDGDq2=fitresult$savefits$sd,
                                            errtype = "tot", fitno=fitresult$savefits$index))

                resultdat$tsnk[4*parindex]       <- tsnk[index]
                resultdat$Nt[4*parindex]         <- Nt[index]
                resultdat$w[4*parindex]          <- datalist[[datanumber]][[2+iset]]$metadata$w
                resultdat$nerr[4*parindex]       <- nerr[index]
                resultdat$iz[4*parindex]         <- iz
                resultdat$icomb[4*parindex]      <- icomb
                resultdat$th[4*parindex]         <- th[index]
                resultdat$errtype[4*parindex]    <- "tot"
                resultdat$DGDq2[4*parindex]      <- fitresult$cdfsummary$t0[1]
                resultdat$dDGDq2[4*parindex]     <- fitresult$cdfsummary$errstat[1]
                resultdat$dat[, 4*parindex]      <- fitresult$cdfbootsamples$bootmean
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

