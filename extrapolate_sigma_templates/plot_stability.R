library(hadron)
library(optparse)

#~ Rscript plot_stability.R --path data/cB211ab.07.64/th6_t48/outputDGammaDq2 --plotpath data/cB211ab.07.64/th6_t48/outputDGammaDq2
#~ Rscript plot_stability.R --path data/cB211ab.07.64/th6_t52/outputDGammaDq2 --plotpath data/cB211ab.07.64/th6_t52/outputDGammaDq2
#~ Rscript plot_stability.R --path data/cB211ab.07.64/th6_t56/outputDGammaDq2 --plotpath data/cB211ab.07.64/th6_t56/outputDGammaDq2

#~ Rscript plot_stability.R --path data/cB211ab.07.64/cov_diag/t48/outputDGammaDq2 --plotpath data/cB211ab.07.64/cov_diag/t48/
#~ Rscript plot_stability.R --path data/cB211ab.07.64/cov_diag/t52/outputDGammaDq2 --plotpath data/cB211ab.07.64/cov_diag/t52/
#~ Rscript plot_stability.R --path data/cB211ab.07.64/cov_diag/t56/outputDGammaDq2 --plotpath data/cB211ab.07.64/cov_diag/t56/
#~ Rscript plot_stability.R --path data/cB211ab.07.64/th6_t56_44/outputDGammaDq2 --plotpath data/cB211ab.07.64/th6_t56_44/

if (TRUE) {
    # set option list
option_list <- list(
make_option(c("--path"), type="character", default="./",
            help = "path to where the data are stored"),
make_option(c("--plotpath"), type="character", default="./",
            help = "path to where the plots are stored"),
make_option(c("--neps"), type="integer", default=10,
            help = "number of epsilons that are considered"),
make_option(c("--nset"), type="integer", default=1,
            help = "number of sets that are considered"),
make_option(c("--nnorm"), type="integer", default=2,
            help = "number of norms that are considered"),
make_option(c("--comment"), type="character", default="",
            help = "comment to put on end of files"),
make_option(c("--calcZ"), action = "store_true", default = FALSE, 
            help = "if given, plot the Z-correlators"),
make_option(c("--inputfile"), type="character", default="",
            help = "inputfile to read in norms")
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
}


givesortednorms <- function(inputfile) {
    input <- read.table(inputfile, fill=TRUE, header=FALSE, comment.char = "[", colClasses=c("character", "character", NULL, NULL))
    nnorms <- as.numeric(input$V2[input$V1=="nnorms"])
    norms <- data.frame(alpha=rep(NA, nnorms), power=rep(NA, nnorms))
    for (i in seq(1, nnorms)) {
        print(input$V2[input$V1==paste0("norm", i, "_am_H*alpha")])
        print(input$V2[input$V1==paste0("norm", i, "_power")])
        norms$alpha[i] <- as.numeric(input$V2[input$V1==paste0("norm", i, "_am_H*alpha")])
        norms$power[i] <- as.numeric(input$V2[input$V1==paste0("norm", i, "_power")])
    }
    norms$legend <- paste("alpha = ", norms$alpha, "pow = ", norms$power)
    return (norms[order(-norms$alpha), ])
}

if (file.exists(opt$inputfile)) {
    legendtext <- givesortednorms(opt$inputfile)$legend
} else {
    legendtext <- seq(1, opt$nnorm)
}


pdf(sprintf("%s/stability_lambda_%s.pdf", opt$plotpath, opt$comment), title="")

for (iset in seq(0, (opt$nset-1))) {
    for (icomb in seq(0, 4)) {
        for (iz in seq(0, 2)) {
            for (ieps in seq(0, (opt$neps-1))) {
                filename <- sprintf("%s/DG_%d_iset_%d_ieps_%d_icomb_%d.dat", opt$path, iz, iset, ieps, icomb)
                title <- sprintf("iz %d iset %d ieps %d icomb %d", iz, iset, ieps, icomb)

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
                    print(title)
                    data <- data[-seq(1, 3*opt$nnorm), ]
                    data$ik <- as.integer(data$ik)
                    data <- data[data$spectreflag == 1, ]

                    try(plot(NA, xlim=c(min(data$lambdalambda_start), max(data$lambdalambda_start)),
                        ylim=c(min(data$rho - data$drho_stat), max(data$rho + data$drho_stat)),
                        xlab="lambda/lambda_start", ylab="rho",
                        main=title,
                        log="x"))

                    for (inorm in seq(1, (opt$nnorm))) {
                        if (inorm == 1) {
                            try(xrange <- c(min(data$lambdalambda_start), max(data$lambdalambda_start)) * c(0.5, 2))
                            try(polygon(x=c(xrange, rev(xrange)),
                                y=c(rep(data$rho[data$ik == inorm & data$resflag == 1] - data$drho_stat[data$ik == inorm & data$resflag == 1], 2),
                                    rep(data$rho[data$ik == inorm & data$resflag == 1] + data$drho_stat[data$ik == inorm & data$resflag == 1], 2)),
                                col=adjustcolor("red", alpha.f=0.2), border=NA))
    
                            try(lines(x=xrange,
                                y=rep(data$rho[data$ik == inorm & data$resflag == 1], 2), col=inorm))
                        }

                        try(plotwitherror(x=data$lambdalambda_start[data$ik == inorm],
                                    y=data$rho[data$ik == inorm],
                                    dy=data$drho_stat[data$ik == inorm],
                                    pch=inorm, col=inorm, rep=TRUE))
                    }
                    legend(x="bottomleft", col=seq(1, opt$nnorm),
                        pch=seq(1, opt$nnorm), legend=legendtext, title="norm")
                }
            }
        }
    }
}

pdf(sprintf("%s/stability_AA0_%s.pdf", opt$plotpath, opt$comment), title="")

for (iset in seq(0, (opt$nset-1))) {
    for (icomb in seq(0, 4)) {
        for (iz in seq(0, 2)) {
            for (ieps in seq(0, (opt$neps-1))) {
                filename <- sprintf("%s/DG_%d_iset_%d_ieps_%d_icomb_%d.dat", opt$path, iz, iset, ieps, icomb)
                title <- sprintf("iz %d iset %d ieps %d icomb %d", iz, iset, ieps, icomb)

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
                    print(title)
                    data <- data[-seq(1, 3*opt$nnorm), ]
                    data$ik <- as.integer(data$ik)
                    data <- data[data$spectreflag == 1, ]

                    try(plot(NA, xlim=c(min(data$AA0_ref), max(data$AA0_ref)),
                        ylim=c(min(data$rho - data$drho_stat), max(data$rho + data$drho_stat)),
                        xlab="[A/A0]_ref", ylab="rho",
                        main=title,
                        log="x"))

                    for (inorm in seq(1, (opt$nnorm))) {
                        if (inorm == 1) {
                            try(xrange <- c(min(data$AA0_ref), max(data$AA0_ref)) * c(0.5, 2))
                            try(polygon(x=c(xrange, rev(xrange)),
                                y=c(rep(data$rho[data$ik == inorm & data$resflag == 1] - data$drho_stat[data$ik == inorm & data$resflag == 1], 2),
                                    rep(data$rho[data$ik == inorm & data$resflag == 1] + data$drho_stat[data$ik == inorm & data$resflag == 1], 2)),
                                col=adjustcolor("red", alpha.f=0.2), border=NA))
    
                            try(lines(x=xrange,
                                y=rep(data$rho[data$ik == inorm & data$resflag == 1], 2), col=inorm))
                        }

                        try(plotwitherror(x=data$AA0_ref[data$ik == inorm],
                                    y=data$rho[data$ik == inorm],
                                    dy=data$drho_stat[data$ik == inorm],
                                    pch=inorm, col=inorm, rep=TRUE))
                    }
                    legend(x="bottomleft", col=seq(1, opt$nnorm),
                        pch=seq(1, opt$nnorm), legend=legendtext, title="norm")
                }
            }
        }
    }
}

pdf(sprintf("%s/reconstruct_kernel%s.pdf", opt$plotpath, opt$comment), title="")

for (iset in seq(0, (opt$nset-1))) {
    for (icomb in seq(0, 4)) {
        for (iz in seq(0, 2)) {
            for (ieps in seq(0, (opt$neps-1))) {
                filename <- sprintf("%s/DG_%d_iset_%d_ieps_%d_icomb_%d.dat", opt$path, iz, iset, ieps, icomb)
                title <- sprintf("iz %d iset %d ieps %d icomb %d", iz, iset, ieps, icomb)

## we use fill=TRUE to ensure that all files are read even if they do not have the same length
## The first lines are details about the norms, three per norm. We delete them.
## spectreflag tells us if we are looking at the stability plot or reconstructed kernel.
## columns 8-18 were only needed for the stability plot, so we do not give them names
                data <- try(read.table(filename, fill=TRUE,
                        col.names=c("ik", "spectreflag", "omega",
                         "kernel", "kernelbar", "kernelbarmkernel", "aM",
                         paste0("column", 8:18))))
                if(!inherits(data, "try-error")) {
                    data <- data[-seq(1, 3*opt$nnorm), ]
                    data$ik <- as.integer(data$ik)
                    data <- data[data$spectreflag == 0, ]


                    for (inorm in seq(1, (opt$nnorm))) {
                        plot(NA, xlim=c(min(data$omega), max(data$omega)),
                            ylim=c(min(c(0, data$kernelbar)), max(data$kernelbar)),
                            xlab="omega", ylab="kernel",
                            main=paste(title, legendtext[inorm]))

                        lines(x=data$omega[data$ik == inorm],
                                y=data$kernel[data$ik == inorm],
                                col=1)

                        lines(x=data$omega[data$ik == inorm],
                                y=data$kernelbar[data$ik == inorm],
                                col=2)
                                
                        legend(x="bottomleft", col=seq(1, 2),
                                pch=c(1, 1), legend=c("kernel", "reconstructed"))
                    }
                }
            }
        }
    }
}


if (opt$calcZ) {
    pdf(sprintf("%s/Z.pdf", opt$plotpath), title="")
    filename <- sprintf("%s/Z.dat", opt$path)
    Z <- try(read.table(filename, fill=TRUE,
                    col.names=c("m", "dm", "iset", "w", "iz", "t", "z", "dz",
                    paste0("column", 9:16))))
                    
    plotwitherror(x=Z$t[Z$iz==0], y=Z$z[Z$iz==0], dy=Z$dz[Z$iz==0],
                    log="y", col=1, pch=1,
                    ylim=c(min(Z$z[Z$z > 0]), max(Z$z)), 
                    xlab="t/a", ylab="Z")
    plotwitherror(x=Z$t[Z$iz==1], y=Z$z[Z$iz==1], dy=Z$dz[Z$iz==1],
                    rep=TRUE, col=2, pch=2)
    plotwitherror(x=Z$t[Z$iz==2], y=Z$z[Z$iz==2], dy=Z$dz[Z$iz==2],
                    rep=TRUE, col=4, pch=4)
    legend(x="topright", legend=c("Z_0", "Z_1", "Z_2"), col=c(1, 2, 4), pch=c(1, 2, 4))
}
