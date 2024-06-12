library(hadron)

givesortednorms <- function(inputfile) {
    input <- read.table(inputfile, fill=TRUE, header=FALSE, comment.char = "[", colClasses=c("character", "character", NULL, NULL))
    nnorms <- as.numeric(input$V2[input$V1=="nnorms"])
    norms <- data.frame(alpha=rep(NA, nnorms), power=rep(NA, nnorms))
    for (i in seq(1, nnorms)) {
#~         print(input$V2[input$V1==paste0("norm", i, "_am_H*alpha")])
#~         print(input$V2[input$V1==paste0("norm", i, "_power")])
        norms$alpha[i] <- as.numeric(input$V2[input$V1==paste0("norm", i, "_am_H*alpha")])
        norms$power[i] <- as.numeric(input$V2[input$V1==paste0("norm", i, "_power")])
    }
    norms$legend <- paste("alpha = ", norms$alpha, "pow = ", norms$power)
    return (norms[order(-norms$alpha), ])
}
makelegendtext <- function(inputfile, nnorm=3){
    if (file.exists(inputfile)) {
        legendtext <- givesortednorms(inputfile)$legend
    } else {
        legendtext <- seq(1, nnorm)
    }
    return(legendtext)
}


#~ pdf(sprintf("%s/stability_lambda_%s.pdf", opt$plotpath, opt$comment), title="")

plot_stability_lambda <- function(path, nset=1, neps=17, nnorm=3, inputfile="", epslist=-1, zlist=-1, comblist=-1, mode="DG") {
    epsseq <- epslist
    if (epslist[1]==-1) epsseq <- seq(0, neps-1)
    combseq <- comblist
    if (comblist[1]==-1) combseq <- seq(0, 4)
    zseq <- zlist
    if (zlist[1]==-1) zseq <- seq(0, 2)
    if (zlist[1]==-1 && mode=="DM") zseq <- seq(0, 3)
    for (iset in seq(0, (nset-1))) {
        for (icomb in combseq) {
            for (iz in zseq) {
                for (ieps in epsseq) {
                    filename <- sprintf("%s/%s_%d_iset_%d_ieps_%d_icomb_%d.dat", path, mode, iz, iset, ieps, icomb)
                    title <- sprintf("iz %d iset %d ieps %d icomb %d", iz, iset, ieps, icomb)
                    legendtext <- makelegendtext(inputfile, nnorm)
    
    ## we use fill=TRUE to ensure that all files are read even if they do not have the same length
    ## The first lines are details about the norms, three per norm. We delete them.
    ## spectreflag tells us if we are looking at the stability plot or reconstructed kernel.
    ## resflag tells us if we are looking at the intermediate steps or the result.
    ## plot the result as a band first.
                    data <- try(read.table(filename, fill=TRUE,
                            col.names=c("ik", "spectreflag", "lambda_start",
                             "lambdalambda_start", "Bnorm", "A0", "AA0_min",
                             "AA0_ref", "AA0", "BnormB_ref", "BnormB",
                             "C_ref", "C", "rho", "drho_stat", "drho_syst",
                             "drho_tot", "resflag")))
                    if(!inherits(data, "try-error")) {
                        data <- data[-seq(1, 3*nnorm), ]
                        data$ik <- as.integer(data$ik)
                        data <- data[data$spectreflag == 1, ]
    
                        try(plot(NA, xlim=c(min(data$lambdalambda_start), max(data$lambdalambda_start)),
                            ylim=c(min(data$rho - data$drho_stat), max(data$rho + data$drho_stat)),
                            xlab="lambda/lambda_start", ylab="rho",
                            main=title,
                            log="x"))
                        
                        legend(x="bottomleft", col=seq(1, nnorm),
                            pch=seq(1, nnorm), legend=legendtext, title="norm")
    
                        for (inorm in seq(0, (nnorm-1))) {
                            if (inorm == 0) {
                                try(xrange <- c(min(data$lambdalambda_start), max(data$lambdalambda_start)) * c(0.5, 2))
                                try(polygon(x=c(xrange, rev(xrange)),
                                    y=c(rep(data$rho[data$ik == inorm & data$resflag == 1] - data$drho_stat[data$ik == inorm & data$resflag == 1], 2),
                                        rep(data$rho[data$ik == inorm & data$resflag == 1] + data$drho_stat[data$ik == inorm & data$resflag == 1], 2)),
                                    col=adjustcolor("red", alpha.f=0.2), border=NA))
        
                                try(lines(x=xrange,
                                    y=rep(data$rho[data$ik == inorm & data$resflag == 1], 2), col=inorm+1))
                            }
    
                            try(plotwitherror(x=data$lambdalambda_start[data$ik == inorm],
                                        y=data$rho[data$ik == inorm],
                                        dy=data$drho_stat[data$ik == inorm],
                                        pch=inorm+1, col=inorm+1, rep=TRUE))
                        }
                    }
                }
            }
        }
    }
}

#~ pdf(sprintf("%s/stability_AA0_%s.pdf", opt$plotpath, opt$comment), title="")
plot_stability_AA0 <- function(path, nset=1, neps=17, nnorm=3, inputfile="", epslist=-1, zlist=-1, comblist=-1, mode="DG") {
    epsseq <- epslist
    if (epslist[1]==-1) epsseq <- seq(0, neps-1)
    combseq <- comblist
    if (comblist[1]==-1) combseq <- seq(0, 4)
    zseq <- zlist
    if (zlist[1]==-1) zseq <- seq(0, 2)
    if (zlist[1]==-1 && mode=="DM") zseq <- seq(0, 3)
    for (iset in seq(0, (nset-1))) {
        for (icomb in combseq) {
            for (iz in zseq) {
                for (ieps in epsseq) {
                    filename <- sprintf("%s/%s_%d_iset_%d_ieps_%d_icomb_%d.dat", path, mode, iz, iset, ieps, icomb)
                    title <- sprintf("iz %d iset %d ieps %d icomb %d", iz, iset, ieps, icomb)
                    legendtext <- makelegendtext(inputfile, nnorm)
    
    ## we use fill=TRUE to ensure that all files are read even if they do not have the same length
    ## The first lines are details about the norms, three per norm. We delete them.
    ## spectreflag tells us if we are looking at the stability plot or reconstructed kernel.
    ## resflag tells us if we are looking at the intermediate steps or the result.
    ## plot the result as a band first.
                    data <- try(read.table(filename, fill=TRUE,
                            col.names=c("ik", "spectreflag", "lambda_start",
                             "lambdalambda_start", "Bnorm", "A0", "AA0_min",
                             "AA0_ref", "AA0", "BnormB_ref", "BnormB",
                             "C_ref", "C", "rho", "drho_stat", "drho_syst",
                             "drho_tot", "resflag")))
                    if(!inherits(data, "try-error")) {
                        data <- data[-seq(1, 3*nnorm), ]
                        data$ik <- as.integer(data$ik)
                        data <- data[data$spectreflag == 1, ]
#~                         print(data)
    
                        try(plot(NA, xlim=c(min(data$AA0_ref), max(data$AA0_ref)),
                            ylim=c(min(data$rho - data$drho_stat), max(data$rho + data$drho_stat)),
                            xlab="[A/A0]_ref", ylab="rho",
                            main=title,
                            log="x"))
                        
                        legend(x="bottomleft", col=seq(1, nnorm),
                            pch=seq(1, nnorm), legend=legendtext, title="norm")
    
                        for (inorm in seq(0, (nnorm-1))) {
                            if (inorm == 0) {
                                try(xrange <- c(min(data$AA0_ref), max(data$AA0_ref)) * c(0.5, 2))
                                try(polygon(x=c(xrange, rev(xrange)),
                                    y=c(rep(data$rho[data$ik == inorm & data$resflag == 1] - data$drho_stat[data$ik == inorm & data$resflag == 1], 2),
                                        rep(data$rho[data$ik == inorm & data$resflag == 1] + data$drho_stat[data$ik == inorm & data$resflag == 1], 2)),
                                    col=adjustcolor("red", alpha.f=0.2), border=NA))
        
                                try(lines(x=xrange,
                                    y=rep(data$rho[data$ik == inorm & data$resflag == 1], 2), col=inorm+1))
                            }
    
                            try(plotwitherror(x=data$AA0_ref[data$ik == inorm],
                                        y=data$rho[data$ik == inorm],
                                        dy=data$drho_stat[data$ik == inorm],
                                        pch=inorm+1, col=inorm+1, rep=TRUE))
                        }
                    }
                }
            }
        }
    }
}

#~ pdf(sprintf("%s/reconstruct_kernel%s.pdf", opt$plotpath, opt$comment), title="")
## with plotnorm we can select a specific norm to be plotted
## -1 means plotting all norms
plot_reconstruction <- function(path, nset=1, neps=17, nnorm=3, plotnorm=c(-1), inputfile="", epslist=-1, zlist=-1, comblist=-1, combine=TRUE, passedtitle="", mode="DG") {
    stopifnot(max(plotnorm) < nnorm)
    if (plotnorm[1] == -1) plotnorm <- seq(1, (nnorm))
    epsseq <- epslist
    if (epslist[1]==-1) epsseq <- seq(0, neps-1)
    combseq <- comblist
    if (comblist[1]==-1) combseq <- seq(0, 4)
    zseq <- zlist
    if (zlist[1]==-1) zseq <- seq(0, 2)
    if (zlist[1]==-1 && mode=="DM") zseq <- seq(0, 3)
    for (iset in seq(0, (nset-1))) {
        for (icomb in combseq) {
            for (iz in zseq) {
                for (ieps in epsseq) {
                    filename <- sprintf("%s/%s_%d_iset_%d_ieps_%d_icomb_%d.dat", path, mode, iz, iset, ieps, icomb)
                    title <- sprintf("iz %d iset %d ieps %d icomb %d \n%s", iz, iset, ieps, icomb, passedtitle)
                    legendtext <- makelegendtext(inputfile, nnorm)
    
    ## we use fill=TRUE to ensure that all files are read even if they do not have the same length
    ## The first lines are details about the norms, three per norm. We delete them.
    ## spectreflag tells us if we are looking at the stability plot or reconstructed kernel.
    ## columns 8-18 were only needed for the stability plot, so we do not give them names
                    data <- try(read.table(filename, fill=TRUE,
                            col.names=c("ik", "spectreflag", "omega",
                             "kernel", "kernelbar", "kernelbarmkernel", "aM",
                             paste0("column", 8:18))))
                    if(!inherits(data, "try-error")) {
                        data <- data[-seq(1, 3*nnorm), ]
                        data$ik <- as.integer(data$ik)
                        data <- data[data$spectreflag == 0, ]
    
    
                        if (combine == TRUE ) {
                            plot(NA, xlim=c(min(data$omega), max(data$omega)),
                                ylim=c(min(c(0, data$kernelbar)), max(data$kernelbar)),
                                xlab="omega", ylab="kernel",
                                main=paste(title))
                            }
                        for (inorm in plotnorm) {
                            if (combine == FALSE ) {
                                plot(NA, xlim=c(min(data$omega), max(data$omega)),
                                    ylim=c(min(c(0, data$kernelbar)), max(data$kernelbar)),
                                    xlab="omega", ylab="kernel",
                                    main=paste(title, legendtext[inorm]))
                                }
    
                            lines(x=data$omega[data$ik == inorm-1],
                                    y=data$kernel[data$ik == inorm-1],
                                    col=1)
    
                            lines(x=data$omega[data$ik == inorm-1],
                                    y=data$kernelbar[data$ik == inorm-1],
                                    col=inorm+1)
                                    
                             if (combine == FALSE) {
                                 legend(x="topright", col=seq(1, inorm+1),
                                    pch=c(1, 1), legend=c("kernel", "reconstructed"))
                                }
                        }
                        if (combine ==TRUE) {
                            legend(x="topright", col=seq(1, nnorm+1),
                                pch=rep(1, nnorm+1), legend=c("kernel", legendtext), title="norm")
                        }
                    }
                }
            }
        }
    }
}


plotZ <- function(path, main="", mode="DG") {
    if(mode!="DG" && mode!="DM") stop("invalid mode, must be DM or DG")
    if(mode=="DG") filename <- sprintf("%s/Z.dat", path)
    if(mode=="DM") filename <- sprintf("%s/EM_Z.dat", path)
    Z <- try(read.table(filename, fill=TRUE,
                    col.names=c("m", "dm", "iset", "w", "iz", "t", "z", "dz",
                    paste0("column", 9:16))))
                    
    plotwitherror(x=Z$t[Z$iz==0], y=Z$z[Z$iz==0], dy=Z$dz[Z$iz==0],
                    log="y", col=1, pch=1,
                    ylim=c(min(Z$z[Z$z > 0]), max(Z$z)), 
                    xlab="t/a", ylab="Z", main=main)
    plotwitherror(x=Z$t[Z$iz==1], y=Z$z[Z$iz==1], dy=Z$dz[Z$iz==1],
                    rep=TRUE, col=2, pch=2)
    plotwitherror(x=Z$t[Z$iz==2], y=Z$z[Z$iz==2], dy=Z$dz[Z$iz==2],
                    rep=TRUE, col=4, pch=4)
    if(mode=="DM") plotwitherror(x=Z$t[Z$iz==3], y=Z$z[Z$iz==3], dy=Z$dz[Z$iz==3],
                    rep=TRUE, col=5, pch=5)
    if(mode=="DG") legend(x="topright", legend=c("Z_0", "Z_1", "Z_2"), col=c(1, 2, 4), pch=c(1, 2, 4))
    if(mode=="DM") legend(x="topright", legend=c("Z_0", "Z_1", "Z_2", "Z_3"), col=c(1, 2, 4, 5), pch=c(1, 2, 4, 5))
}


#~ pdf(sprintf("%s/stability_AA0_%s.pdf", opt$plotpath, opt$comment), title="")
plot_stability_AA0_sets <- function(pathlist, nset=1, neps=17, nnorm=3, inputfile="", epslist=-1, zlist=-1, comblist=-1, mode="DG", comments=c(), normseq=c(1)) {
    mycolours <- c("red", "blue", "green", "cyan", "pink", "#D2691E", "#556B2F", rep("black", 100))
    epsseq <- epslist
    if (epslist[1]==-1) epsseq <- seq(0, neps-1)
    combseq <- comblist
    if (comblist[1]==-1) combseq <- seq(0, 4)
    zseq <- zlist
    if (zlist[1]==-1) zseq <- seq(0, 2)
    if (zlist[1]==-1 && mode=="DM") zseq <- seq(0, 3)
    if(length(normseq)>1) stop("you are plotting more than one norm. This is not possible")
    inorm <- normseq[1]
    if(length(comments)!=length(pathlist)) print("length of comments and pathlist does not match")
    for (iset in seq(0, (nset-1))) {
        for (icomb in combseq) {
            for (iz in zseq) {
                for (ieps in epsseq) {
                    filenamelist <- c()
                    data <- list()
                    do_continue <- T
                    for(ipath in seq_along(pathlist)) {
                        filenamelist[ipath] <- sprintf("%s/%s_%d_iset_%d_ieps_%d_icomb_%d.dat", pathlist[ipath], mode, iz, iset, ieps, icomb)
                        data[[ipath]] <- try(read.table(filenamelist[ipath], fill=TRUE,
                                col.names=c("ik", "spectreflag", "lambda_start",
                                 "lambdalambda_start", "Bnorm", "A0", "AA0_min",
                                 "AA0_ref", "AA0", "BnormB_ref", "BnormB",
                                 "C_ref", "C", "rho", "drho_stat", "drho_syst",
                                 "drho_tot", "resflag")))
                        do_continue <- do_continue && !inherits(data[[ipath]], "try-error")
                    }
                    title <- sprintf("iz %d iset %d ieps %d icomb %d", iz, iset, ieps, icomb)
                    legendtext <- makelegendtext(inputfile, nnorm)
    
    ## we use fill=TRUE to ensure that all files are read even if they do not have the same length
    ## The first lines are details about the norms, three per norm. We delete them.
    ## spectreflag tells us if we are looking at the stability plot or reconstructed kernel.
    ## resflag tells us if we are looking at the intermediate steps or the result.
    ## plot the result as a band first.
                    
                    if(do_continue) {
                        boundx <- c()
                        minimay <- c()
                        maximay <- c()
                        for(ipath in seq_along(pathlist)) {
                            data[[ipath]] <- data[[ipath]][-seq(1, 3*nnorm), ]
                            data[[ipath]]$ik <- as.integer(data[[ipath]]$ik)
                            data[[ipath]] <- data[[ipath]][data[[ipath]]$spectreflag == 1, ]
                            print(ipath)
                            print(range(data[[ipath]]$rho[data[[ipath]]$ik == inorm]))
                            boundx <- append(boundx, data[[ipath]]$AA0_ref[data[[ipath]]$ik == inorm])
                            minimay <- append(minimay, data[[ipath]]$rho[data[[ipath]]$ik == inorm] - data[[ipath]]$drho_stat[data[[ipath]]$ik == inorm])
                            maximay <- append(maximay, data[[ipath]]$rho[data[[ipath]]$ik == inorm] + data[[ipath]]$drho_stat[data[[ipath]]$ik == inorm])
                        }
                        try(plot(NA, xlim=range(boundx),
                            ylim=c(min(minimay), max(maximay)),
                            xlab="[A/A0]_ref", ylab="rho",
                            main=title,
                            log="x"))
                        
                        legend(x="bottomleft",
                            col=mycolours[seq_along(pathlist)],
                            pch=seq_along(pathlist),
                            legend=comments, 
                            title=paste("norm", normseq), border=NA)
    
                        for(ipath in seq_along(pathlist)) {
                            try(xrange <- c(min(data[[ipath]]$AA0_ref), max(data[[ipath]]$AA0_ref)) * c(0.01, 20))
                            try(polygon(x=c(xrange, rev(xrange)),
                                y=c(rep(data[[ipath]]$rho[data[[ipath]]$ik == inorm & data[[ipath]]$resflag == 1] - data[[ipath]]$drho_stat[data[[ipath]]$ik == inorm & data[[ipath]]$resflag == 1], 2),
                                    rep(data[[ipath]]$rho[data[[ipath]]$ik == inorm & data[[ipath]]$resflag == 1] + data[[ipath]]$drho_stat[data[[ipath]]$ik == inorm & data[[ipath]]$resflag == 1], 2)),
                                col=adjustcolor(mycolours[ipath], alpha.f=0.2), border=NA))
                            try(lines(x=xrange,
                                y=rep(data[[ipath]]$rho[data[[ipath]]$ik == inorm & data[[ipath]]$resflag == 1], 2), col=mycolours[ipath]))
                            try(plotwitherror(x=data[[ipath]]$AA0_ref[data[[ipath]]$ik == inorm],
                                        y=data[[ipath]]$rho[data[[ipath]]$ik == inorm],
                                        dy=data[[ipath]]$drho_stat[data[[ipath]]$ik == inorm],
                                        pch=ipath, col=mycolours[ipath], rep=TRUE))
                        }
                    }
                }
            }
        }
    }
}
