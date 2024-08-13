library("hadron")
source("/home/gross/Documents/heavymesons/scripts_christiane/extrapolate_sigma_templates/read_in_binary.R")



## bootsample = bootsample - (mean - bootsample) * (errtotmean / bootstat - 1)
increasedistancebootstrap <- function(sample, mean, ratio) (sample - (mean - sample) * (ratio - 1))


## TODO: different finite volume effects possible
determineDGDq2_contlim <- function(resultpathlist, filenames, tsnk, Nt, th,
        nerr, amin, savename, fitfn, par.guess, isets=c(0), comblist=c(0), zlist=c(0, 1, 2, 3),
        epsuplim=0, epslowlim=0, errors=c("stat"), volumetablelist = c(""),
        doplot=TRUE, NDG = 4, maxz=4, bsamples=1000, afm, mdsgev, numberspacings=4, 
        cols=c("blue", "green", "#D2691E", "#556B2F", "cyan", "pink", "red", rep("black", 100)),
        legendargs=list(), neps=-1) {

    ## all members of errors have to be one of "stat", "sys", "vol", "tot" and no duplicates
    ## for calculating tot we have to calculate sys and vol as well
    for (i in seq(1, length(errors))){
        stopifnot(errors[i] == "stat" || errors[i] == "sys" || errors[i] == "vol" || errors[i] == "tot" )
    }

    stopifnot(!("tot" %in% errors && (!("sys" %in% errors) || !("vol" %in% errors))))

    ## detect some errors
    stopifnot(length(unique(errors)) == length(errors))
    stopifnot(length(resultpathlist) == length(afm))
    stopifnot(length(resultpathlist) == length(mdsgev))
    stopifnot(length(resultpathlist) == numberspacings)
    stopifnot(length(filenames)/numberspacings == length(tsnk))
    stopifnot(length(filenames)/numberspacings == length(Nt))
    stopifnot(length(filenames)/numberspacings == length(th))
    stopifnot(length(filenames)/numberspacings == length(nerr))
    
    ## and some probable mistakes
    if(NDG!=4 && NDG!=5) print("Are you sure about your value of NDG?")
    
    finitevolumeall <- list()
    if ("vol" %in% errors || "tot" %in% errors) {
        stopifnot(length(volumetablelist)==numberspacings)
        for (i in seq(1, numberspacings)) {
            stopifnot(file.exists(volumetablelist[i]))
    ## read in table with finite volume results
            finitevolumeall[[i]] <- read.table(volumetablelist[i], header=TRUE)
        }
    }

    if(neps==-1) neps=epsuplim-epslowlim
    result <- data.frame(tsnk=NA, Nt=NA, w=NA, nerr=NA, iz=NA, icomb=NA, th = NA, epsindex=NA, eps=NA, DGDq2=NA, dDGDq2=NA, errtype=NA)

    resultdat <- list(tsnk=c(), Nt=c(), w=c(), nerr=c(), iz=c(), 
    icomb=c(), th=c(), epsindex=c(), eps=c(), errtype = c(), DGDq2 = c(), dDGDq2 = c(), 
    dat = array(NA, dim=c(bsamples, length(tsnk)*4*length(isets)*length(comblist)*maxz*neps)))
    
    fitlist <- list()
    namefitlist <- seq(1, length(tsnk)*4*length(isets)*length(comblist)*maxz*neps)
    
    legendargs[["col"]]=cols[1:numberspacings]
    legendargs[["pch"]]=rep(1, numberspacings)

    print(seq(epslowlim, epslowlim+neps-1))

    parindex <- 1
    datalist <- list()
for (index in seq(1, length(filenames)/numberspacings)){
    for (datanumber in seq(1, numberspacings)) {
#~         print(filenames[index+(datanumber-1)*length(tsnk)])
        datalist[[datanumber]] <- try(read_in_DGDq2(filename=filenames[index+(datanumber-1)*length(th)], write=FALSE, resultpath=resultpathlist[[datanumber]], NDG = NDG))
    }
for (iset in isets) {
    mymasks <- list(rep(NA, 4), rep(NA, 4), rep(NA, 4), rep(NA, 4), rep(NA, 4))
    ## epsilon are counted in C-style numbering
    for(eps in seq(epslowlim, epslowlim+neps-1)) {
    for (iz in zlist) {
        for(icomb in comblist) {
            do_continue <- TRUE
            for (datanumber in seq(1, numberspacings)) {
                do_continue <- do_continue && !inherits(data, "try-error")
            }
            if(do_continue) {

            DGamma <- c()
            dDGamma <- c()
            dDGammasys <- c()
            bs <- array(NA, dim=c(datalist[[1]][[2+iset]][[3]]$bootnumber[[1]], numberspacings))
            xlist <- list()
            bootsys <- array(NA, dim=c(datalist[[1]][[2+iset]][[3]]$bootnumber[[1]], numberspacings))
            bootvol <- array(NA, dim=c(datalist[[1]][[2+iset]][[3]]$bootnumber[[1]], numberspacings))
            boottot <- array(NA, dim=c(datalist[[1]][[2+iset]][[3]]$bootnumber[[1]], numberspacings))
            

            for (datanumber in seq(1, numberspacings)) {
            ## set number of epsilons: default epslowlim:epsuplim, if neps is given, epslowlim:(epslowlim+neps)
            stopifnot(datalist[[datanumber]][[2+iset]]$metadata$neps >= epslowlim+neps)
            if(neps!=-1) epsuplim <- datalist[[datanumber]][[2+iset]]$metadata$neps - (epslowlim+neps)
            ## determine which part of data to use, extract needed elements
            namesel <- paste0("id", ((iz)*5 + icomb)*datalist[[datanumber]][[2+iset]]$metadata$neps + eps, "ieps", eps, "icomb", icomb, "iz", iz)
            title=paste("iset", iset, "iz", iz, "icomb", icomb, "tsnk", tsnk[index], "Nt", Nt[index], "th", th[index], "eps", eps, "nerr", nerr[index])
            if(datanumber==1) print(title)

            DGamma[datanumber] <- unlist(sapply(X=datalist[[datanumber]][[2+iset]][namesel], FUN=getElement, name="DGDq2mean")[1, ], use.names=F)
            dDGammasys[datanumber] <- unname(unlist(sapply(X=datalist[[datanumber]][[2+iset]][namesel], FUN=getElement, name="sys")[1, ], use.names=F))
            dDGamma[datanumber] <- unlist(sapply(X=datalist[[datanumber]][[2+iset]][namesel], FUN=getElement, name="DGDq2sd")[1, ], use.names=F)
            DGammaboot <- sapply(X=datalist[[datanumber]][[2+iset]][namesel], FUN=getElement, name="DGDq2boots")[1, ]
            
            bs[, datanumber] <- array(unlist(DGammaboot, use.names=F), dim=c(datalist[[datanumber]][[2+iset]][[namesel[1]]]$bootnumber[[1]], 1))
            xlist[[datanumber]] <- datalist[[datanumber]][[2+iset]]$epsilons[(epslowlim+1):(datalist[[datanumber]][[2+iset]]$metadata$neps-epsuplim)]
            
            ## determine min A/A0_ref
            AA0ar <- array(NA, dim=c(datalist[[datanumber]][[2+iset]]$metadata$nnorm, 1))
            for (norm in seq(1, datalist[[datanumber]][[2+iset]]$metadata$nnorm)) {
                A0ABCW_ref <- sapply(X=datalist[[datanumber]][[2+iset]][namesel], FUN=getElement, name="A0ABCW_ref")[1, ]
                A0ABCW_ref_ar <- t(array(unlist(A0ABCW_ref, use.names=F), dim=c(6, 1)))
                AA0ar[norm, ] <- A0ABCW_ref_ar[, 2] / A0ABCW_ref_ar[, 1]
            }
            AA0 <- apply(X=AA0ar, MARGIN=2, FUN=min)

            ## determine masks for fit
            ## TODO
            if(iz != maxz-1) mymasks[[iz+1]][datanumber] <- (AA0 < amin)
            if (iz==maxz-1) {
                if(is.na(mymasks[[0+1]][datanumber])) {
                    mymasks[[maxz]][datanumber] <- 0
                } else {
                    mymasks[[maxz]][datanumber] <- mymasks[[0+1]][datanumber]
                    for (zmask in seq(1, maxz-2)) {
                        mymasks[[maxz]][datanumber] <- mymasks[[maxz]][datanumber] & mymasks[[zmask+1]][datanumber]
                    }
                }
            }
        }

            ## error only from statistics
            if ("stat" %in% errors) {
                ## extrapolate a^2 to zero
                fitresult <- try(bootstrap.nlsfit(x=afm^2, y=DGamma,
                                    bsamples=bs, fn=fitfn, par.guess=par.guess,
                                    mask=mymasks[[iz+1]], verbose=(iz==3)))
                
                if(!inherits(fitresult, "try-error")) {
                if (doplot) try(plot(fitresult, main=paste(title, "only stat"), xlab="a[fm^2]", ylab="DG/Dw2", xlim=c(0, max(afm^2)), plot.range=c(-0.001, 0.007)))
                if (doplot) try(plotwitherror(x=0, y=fitresult$t0[1], dy=fitresult$se[1], rep=TRUE, col="red", cex=1, lwd=3))

                ## save result
                ## save combined result with index -1
                result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=datalist[[datanumber]][[2+iset]]$metadata$w,
                                            nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], epsindex=eps, eps=datalist[[datanumber]][[2+iset]]$epsilons[eps+1], 
                                            DGDq2=fitresult$t0[1], dDGDq2=fitresult$se[1],
                                            errtype = "stat"))
                resultdat$tsnk[4*parindex-3]       <- tsnk[index]
                resultdat$Nt[4*parindex-3]         <- Nt[index]
                resultdat$w[4*parindex-3]          <- datalist[[1]][[2+iset]]$metadata$w
                resultdat$nerr[4*parindex-3]       <- nerr[index]
                resultdat$iz[4*parindex-3]         <- iz
                resultdat$icomb[4*parindex-3]      <- icomb
                resultdat$th[4*parindex-3]         <- th[index]
                resultdat$epsindex[4*parindex-3]   <- eps
                resultdat$eps[4*parindex-3]        <- datalist[[datanumber]][[2+iset]]$epsilons[eps+1]
                resultdat$errtype[4*parindex-3]    <- "stat"
                resultdat$DGDq2[4*parindex-3]      <- fitresult$t0[1]
                resultdat$dDGDq2[4*parindex-3]     <- fitresult$se[1]
                resultdat$dat[, 4*parindex-3]      <- fitresult$t[, 1]
                fitlist[[4*parindex-3]]            <- fitresult
                namefitlist[4*parindex-3]          <- paste(title, "stat")


                } else {
                    result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=datalist[[datanumber]][[2+iset]]$metadata$w, nerr=nerr[index],
                                            iz=iz, icomb=icomb, th=th[index], epsindex=eps, eps=datalist[[datanumber]][[2+iset]]$epsilons[eps+1], DGDq2=NA, dDGDq2=NA, errtype=NA))
                    fitlist[[4*parindex-3]] <- NA
                }
            } else fitlist[[4*parindex-3]] <- NA

#~             ## extrapolate with systematical error from HLT
            if ("sys" %in% errors) {
            for (datanumber in seq(1, numberspacings)) {
                ratio <- sqrt(dDGamma[[datanumber]]^2 + dDGammasys[[datanumber]]^2)/dDGamma[[datanumber]]
                bootsys[, datanumber] <- bs[, datanumber]*(-ratio) + array(rep(DGamma[[datanumber]]*(1+ratio), each=bsamples), dim=c(bsamples, length(ratio)))
                    }
                        
#~ #                print(dDGamma_sys/dDGamma)
                
                ## extrapolate a^2 to zero
                fitresult <- try(bootstrap.nlsfit(x=afm^2, y=DGamma,
                                    bsamples=bootsys, fn=fitfn, par.guess=par.guess,
                                    mask=mymasks[[iz+1]], verbose=(iz==3)))
                
                if(!inherits(fitresult, "try-error")) {
                if (doplot) try(plot(fitresult, main=paste(title, "stat + sys"), xlab="a[fm^2]", ylab="DG/Dw2", xlim=c(0, max(afm^2)), plot.range=c(-0.001, 0.007)))
                if (doplot) try(plotwitherror(x=0, y=fitresult$t0[1], dy=fitresult$se[1], rep=TRUE, col="red", cex=1, lwd=3))

                ## save result
                ## save combined result with index -1
                result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=datalist[[datanumber]][[2+iset]]$metadata$w,
                                            nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], epsindex=eps, eps=datalist[[datanumber]][[2+iset]]$epsilons[eps+1], 
                                            DGDq2=fitresult$t0[1], dDGDq2=fitresult$se[1],
                                            errtype = "sys"))
                resultdat$tsnk[4*parindex-2]       <- tsnk[index]
                resultdat$Nt[4*parindex-2]         <- Nt[index]
                resultdat$w[4*parindex-2]          <- datalist[[1]][[2+iset]]$metadata$w
                resultdat$nerr[4*parindex-2]       <- nerr[index]
                resultdat$iz[4*parindex-2]         <- iz
                resultdat$icomb[4*parindex-2]      <- icomb
                resultdat$th[4*parindex-2]         <- th[index]
                resultdat$epsindex[4*parindex-2]   <- eps
                resultdat$eps[4*parindex-2]        <- datalist[[datanumber]][[2+iset]]$epsilons[eps+1]
                resultdat$errtype[4*parindex-2]    <- "sys"
                resultdat$DGDq2[4*parindex-2]      <- fitresult$t0[1]
                resultdat$dDGDq2[4*parindex-2]     <- fitresult$se[1]
                resultdat$dat[, 4*parindex-2]      <- fitresult$t[, 1]
                fitlist[[4*parindex-2]]            <- fitresult
                namefitlist[4*parindex-2]          <- paste(title, "sys")

                } else {
                    result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=datalist[[datanumber]][[2+iset]]$metadata$w, nerr=nerr[index],
                                            iz=iz, icomb=icomb, th=th[index], epsindex=eps, eps=datalist[[datanumber]][[2+iset]]$epsilons[eps+1], DGDq2=NA, dDGDq2=NA, errtype=NA))
                    fitlist[[4*parindex-2]] <- NA
                }
            } else fitlist[[4*parindex-2]] <- NA

            ## extrapolate with finite volume error
            if ("vol" %in% errors) {
                for (datanumber in seq(1, numberspacings)) {
                    volume <- finitevolumeall[[datanumber]][finitevolumeall[[datanumber]]$th == th[index] & finitevolumeall[[datanumber]]$iz == iz, ]
                    stopifnot(datalist[[datanumber]][[2+iset]]$epsilons == volume$epsilon)
                    
    
                ratio <- sqrt(dDGamma[[datanumber]]^2 + volume$Delta[eps+1]^2)/dDGamma[[datanumber]]
                bootvol[, datanumber] <- bs[, datanumber]*(-ratio) + array(rep(DGamma[[datanumber]]*(1+ratio), each=bsamples), dim=c(bsamples, length(ratio)))
                }
                    
                fitresult <- try(bootstrap.nlsfit(x=afm^2, y=DGamma,
                                    bsamples=bootvol, fn=fitfn, par.guess=par.guess,
                                    mask=mymasks[[iz+1]], verbose=(iz==3)))
                
                if(!inherits(fitresult, "try-error")) {
                if (doplot) try(plot(fitresult, main=paste(title, "stat + vol"), xlab="a[fm^2]", ylab="DG/Dw2", xlim=c(0, max(afm^2)), plot.range=c(-0.001, 0.007)))
                if (doplot) try(plotwitherror(x=0, y=fitresult$t0[1], dy=fitresult$se[1], rep=TRUE, col="red", cex=1, lwd=3))

                ## save result
                ## save combined result with index -1
                result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=datalist[[datanumber]][[2+iset]]$metadata$w,
                                            nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], epsindex=eps, eps=datalist[[datanumber]][[2+iset]]$epsilons[eps+1], 
                                            DGDq2=fitresult$t0[1], dDGDq2=fitresult$se[1],
                                            errtype = "vol"))
                resultdat$tsnk[4*parindex-1]       <- tsnk[index]
                resultdat$Nt[4*parindex-1]         <- Nt[index]
                resultdat$w[4*parindex-1]          <- datalist[[1]][[2+iset]]$metadata$w
                resultdat$nerr[4*parindex-1]       <- nerr[index]
                resultdat$iz[4*parindex-1]         <- iz
                resultdat$icomb[4*parindex-1]      <- icomb
                resultdat$th[4*parindex-1]         <- th[index]
                resultdat$epsindex[4*parindex-1]   <- eps
                resultdat$eps[4*parindex-1]        <- datalist[[datanumber]][[2+iset]]$epsilons[eps+1]
                resultdat$errtype[4*parindex-1]    <- "vol"
                resultdat$DGDq2[4*parindex-1]      <- fitresult$t0[1]
                resultdat$dDGDq2[4*parindex-1]     <- fitresult$se[1]
                resultdat$dat[, 4*parindex-1]      <- fitresult$t[, 1]
                fitlist[[4*parindex-1]]            <- fitresult
                namefitlist[4*parindex-1]          <- paste(title, "vol")


                } else {
                    result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=datalist[[datanumber]][[2+iset]]$metadata$w, nerr=nerr[index],
                                            iz=iz, icomb=icomb, th=th[index], epsindex=eps, eps=datalist[[datanumber]][[2+iset]]$epsilons[eps+1], DGDq2=NA, dDGDq2=NA, errtype=NA))
                    fitlist[[4*parindex-1]] <- NA
                }
            } else fitlist[[4*parindex-1]] <- NA

            ## extrapolate with finite volume and systematic error
            if ("tot" %in% errors) {
                for (datanumber in seq(1, numberspacings)) {
                    volume <- finitevolumeall[[datanumber]][finitevolumeall[[datanumber]]$th == th[index] & finitevolumeall[[datanumber]]$iz == iz, ]
                    stopifnot(datalist[[datanumber]][[2+iset]]$epsilons == volume$epsilon)
                    
                    ratio <- sqrt(dDGamma[[datanumber]]^2 + dDGammasys[[datanumber]]^2 + volume$Delta[eps+1]^2)/dDGamma[[datanumber]]
                    boottot[, datanumber] <- bs[, datanumber]*(-ratio) + array(rep(DGamma[[datanumber]]*(1+ratio), each=bsamples), dim=c(bsamples, length(ratio)))
                }
                
                
                fitresult <- try(bootstrap.nlsfit(x=afm^2, y=DGamma,
                                    bsamples=boottot, fn=fitfn, par.guess=par.guess,
                                    mask=mymasks[[iz+1]], verbose=(iz==3)))
                
                if(!inherits(fitresult, "try-error")) {
                if (doplot) try(plot(fitresult, main=paste(title, "stat + sys + vol"), xlab="a[fm^2]", ylab="DG/Dw2", xlim=c(0, max(afm^2)), plot.range=c(-0.001, 0.007)))
                if (doplot) try(plotwitherror(x=0, y=fitresult$t0[1], dy=fitresult$se[1], rep=TRUE, col="red", cex=1, lwd=3))

                ## save result
                ## save combined result with index -1
                result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=datalist[[datanumber]][[2+iset]]$metadata$w,
                                            nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], epsindex=eps, eps=datalist[[datanumber]][[2+iset]]$epsilons[eps+1], 
                                            DGDq2=fitresult$t0[1], dDGDq2=fitresult$se[1],
                                            errtype = "vol"))
                resultdat$tsnk[4*parindex]       <- tsnk[index]
                resultdat$Nt[4*parindex]         <- Nt[index]
                resultdat$w[4*parindex]          <- datalist[[1]][[2+iset]]$metadata$w
                resultdat$nerr[4*parindex]       <- nerr[index]
                resultdat$iz[4*parindex]         <- iz
                resultdat$icomb[4*parindex]      <- icomb
                resultdat$th[4*parindex]         <- th[index]
                resultdat$epsindex[4*parindex]   <- eps
                resultdat$eps[4*parindex]        <- datalist[[datanumber]][[2+iset]]$epsilons[eps+1]
                resultdat$errtype[4*parindex]    <- "tot"
                resultdat$DGDq2[4*parindex]      <- fitresult$t0[1]
                resultdat$dDGDq2[4*parindex]     <- fitresult$se[1]
                resultdat$dat[, 4*parindex]      <- fitresult$t[, 1]
                fitlist[[4*parindex]]            <- fitresult
                namefitlist[4*parindex]          <- paste(title, "tot")



                } else {
                    result <- rbind(result, data.frame(tsnk=tsnk[index], Nt=Nt[index], w=datalist[[datanumber]][[2+iset]]$metadata$w,
                                                nerr=nerr[index], iz=iz, icomb=icomb, th=th[index], epsindex=eps, eps=datalist[[datanumber]][[2+iset]]$epsilons[eps+1], DGDq2=NA, dDGDq2=NA,
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
                resultdat$epsindex[4*parindex]   <- NA
                resultdat$eps[4*parindex]        <- NA
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
                resultdat$epsindex[4*parindex]   <- NA
                resultdat$eps[4*parindex]        <- NA
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
}
result <- result[-1, ]
result <- result[ ! (result$iz==0 & (result$icomb==3 | result$icomb==4)), ]
print(warnings())
write.table(result, paste0(savename, ".csv"), row.names=F, col.names=T)
saveRDS(resultdat, paste0(savename, ".RDS"))
names(fitlist) <- namefitlist
return(invisible(fitlist))
}

