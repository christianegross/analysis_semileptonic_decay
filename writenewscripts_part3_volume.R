parent <- "~/Documents/heavymesons/data/newinput"
folders <- c("cB211.07.64", "cB211.07.96", "cC211.06.80", "cD211.054.96", "cE211.044.112", "cB211.07.48_300", "cB211.07.48_400", "cB211.07.64_48_36")
nameshort <- c("B64", "B96", "C80", "D96", "E112", "B48_300", "B48_400", "B64_48_36")
subfolders  <- c("th1", "th2", "th3", "th4", "th5", "th6", "th7", "th8", "th9", "th9.5")

infos <- read.table("parameters_input_files.csv", row.names=1, header=TRUE, sep=",", comment.char = "#")
#~ print(infos)
names(infos)


## TODO: write volume effects



fileDG <- sprintf("%s/%s-volume-DGDq2.Rmd", parent, 4*length(folders)+4)
fileDM <- sprintf("%s/%s-volume-DMDq2.Rmd", parent, 4*length(folders)+4)

## volume DG

cat("```{r setupsummaryvolume, include=FALSE}", 
"knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, results = FALSE)", 
"source(\"/home/gross/Documents/heavymesons/scripts_christiane/extrapolate_sigma_templates/functions_stability_plots.R\")", 
"source(\"/home/gross/Documents/heavymesons/scripts_christiane/extrapolate_sigma_templates/calc_DGDq2.R\")", 
"library(\"hadron\")", 
"zlist <- c(3, 0, 1, 2)", 
"errlist <- c(\"stat\")", 
"sumup <-FALSE", 
"```",
"", 
"# Differential Decay Rates - Volume",
file=fileDG, append=FALSE, sep="\n")


for(kernel in c("sigmoid", "erf")){
    for(channel in c("cd", "cs")) {
        cat(sprintf("## %s %s\n", channel, kernel), 
sprintf("```{r decayrate_%s_volume_%s}", channel, kernel), 
sprintf("B64_56_44 <- read.table(\"tables/DG_complete_cB64_converted_%s_finvol_%s.csv\")", channel, kernel), 
sprintf("B64_48_36 <- read.table(\"tables/DG_complete_cB64_48_36_converted_%s_finvol_%s.csv\")", channel, kernel), 
sprintf("B96 <- read.table(\"tables/DG_complete_cB96_converted_%s_finvol_%s.csv\")", channel, kernel), 
sprintf("B48 <- read.table(\"tables/DG_complete_cB48_400_converted_%s_finvol_%s.csv\")", channel, kernel), 
"", 
"", 
sprintf("B64_56_44dat <- readRDS(\"tables/DG_complete_cB64_converted_%s_finvol_%s.RDS\")", channel, kernel), 
sprintf("B64_48_36dat <- readRDS(\"tables/DG_complete_cB64_48_36_converted_%s_finvol_%s.RDS\")", channel, kernel), 
sprintf("B96dat <- readRDS(\"tables/DG_complete_cB96_converted_%s_finvol_%s.RDS\")", channel, kernel), 
sprintf("B48dat <- readRDS(\"tables/DG_complete_cB48_400_converted_%s_finvol_%s.RDS\")", channel, kernel), 
"", 
"B64_56_44 <- na.omit(B64_56_44)", 
"B64_48_36 <- na.omit(B64_48_36)", 
"B96 <- na.omit(B96)", 
"B48 <- na.omit(B48)", 
"", 
"B64_56_44 <- B64_56_44[order(B64_56_44$iz, B64_56_44$errtype, B64_56_44$icomb, B64_56_44$q), ]", 
"B64_48_36 <- B64_48_36[order(B64_48_36$iz, B64_48_36$errtype, B64_48_36$icomb, B64_48_36$q), ]", 
"B96 <- B96[order(B96$iz, B96$errtype, B96$icomb, B96$q), ]", 
"B48 <- B48[order(B48$iz, B48$errtype, B48$icomb, B48$q), ]", 
"", 
"B64_56_44$L <- rep(64, length(B64_56_44$iz))", 
"B64_48_36$L <- rep(64, length(B64_48_36$iz))", 
"B96$L <- rep(96, length(B96$iz))", 
"B48$L <- rep(48, length(B48$iz))", 
"", 
"", 
"tablelist <- list(B48, B64_48_36, B64_56_44, B96)", 
"enslist <- list(B48dat, B64_48_36dat, B64_56_44dat, B96dat)", 
"", 
"", 
"", 
"fnlin <- function(par, x, boot.R, ...) par[1] + par[2] * x", 
"fncon <- function(par, x, boot.R, ...) par[1]", 
"", 
"fitfn <- fnlin", 
"par.guess <- c(1, 1)", 
"comment <- \"linear\"", 
"", 
"res <- data.frame(th=c(), errtype=c(), iz=c(), lim=c(), dlim=c(), sloipe=c(), dslope=c(), chi=c(), q=c(), ",
"                  ratio64=c(), dratio64=c(), ratio5644=c(), dratio5644=c(), ratio4836=c(), dratio4836=c())", 
" ## ratios: (48, 36) / (56, 44); 48/64; 64/96",
" ## trend=ratio larger/smaller 1",
"", 
"", 
"pcol <- col2rgb(\"gray\", alpha=TRUE)/255 ", 
"pcol[4] <- 0.65", 
"pcol <- rgb(red=pcol[1],green=pcol[2],blue=pcol[3],alpha=pcol[4])", 
"", 
"confcounter <- 1", 
"", 
"", 
"", 
"for(th in c(seq(1, 9), 9.5)) {", 
"  for(iz in zlist) {", 
"    for (err in errlist){", 
"      ", 
"      title <- sprintf(\"th%s q^2 %.3f GeV^2 iz %d errtype %s\", th, B64_56_44$q[abs(B64_56_44$th - th) < 1e-2][1]^2, iz, err)", 
"      print(title)", 
"      bsamples <- array(NA, dim=c(1000, length(enslist)))", 
"      y <- c()", 
"      dy <- c()", 
"      x <- c()", 
"      for (index in seq(1, length(enslist))) {", 
"        if(!is.na(tablelist[[index]]$q[abs(tablelist[[index]]$th - th) < 1e-2][1]^2)) {", 
"          mask <- abs(enslist[[index]]$th - th) < 1e-2 & enslist[[index]]$iz == iz & enslist[[index]]$errtype==err", 
"          masktable <- abs(tablelist[[index]]$th - th) < 1e-2 & tablelist[[index]]$iz == iz & tablelist[[index]]$errtype==err", 
"          mask[is.na(mask)] <- F", 
"          y[index] <- enslist[[index]]$DGDq2gev[mask]", 
"          dy[index] <- enslist[[index]]$dDGDq2gev[mask]", 
"          try(bsamples[, index] <- enslist[[index]]$datgev[, mask])", 
"          x[index] <- tablelist[[index]]$L[1]", 
"        }", 
"      }", 
"      ", 
"      maskB64_56_44 <- B64_56_44$iz==iz & B64_56_44$icomb==0 & B64_56_44$errtype==err & B64_56_44$th==th", 
"      maskB64_48_36 <- B64_48_36$iz==iz & B64_48_36$icomb==0 & B64_48_36$errtype==err & B64_48_36$th==th", 
"      maskB96 <- B96$iz==iz & B96$icomb==0 & B96$errtype==err & B96$th==th", 
"      maskB48 <- B48$iz==iz & B48$icomb==0 & B48$errtype==err & B48$th==th", 
"      ", 
"      fitresult <- try(bootstrap.nlsfit(x=1/x, y=y, bsamples=bsamples, fn=fitfn, par.guess=par.guess))", 
"      ", 
"      if(!inherits(x = fitresult, what=\"try-error\")){", 
"        plot(fitresult,", 
sprintf("             main=paste(\"%s %s diff. decay rate , \", err, \", z =\", iz, \"th\", th), ", channel, kernel), 
"             xlab=\"1/L\", ylab=\"24 pi^3 DGamma / Dq^2 [GeV^-3]\", xaxt=\"n\")", 
"        plotwitherror(x=1/B64_56_44$L[maskB64_56_44], y=B64_56_44$DGDq2gev[maskB64_56_44], dy=B64_56_44$dDGDq2gev[maskB64_56_44], rep=T, pch=1)", 
"        ratio64 <- bsamples[, 2] / bsamples[, 3]",
"        ratio5644 <- bsamples[, 3] / bsamples[, 4]",
"        ratio4836 <- bsamples[, 1] / bsamples[, 2]",
"        newline <- data.frame(th=th, errtype=err, iz=iz, ", 
"                              lim=fitresult$t0[1], dlim=fitresult$se[1], ", 
"                              slope=fitresult$t0[2], dslope=fitresult$se[2], ", 
"                              chi=fitresult$chisqr / fitresult$dof, q=B64_56_44$q[abs(B64_56_44$th - th) < 1e-2][1],",
"                              ratio64 = mean(ratio64), dratio64=sd(ratio64),", 
"                              ratio5644 = mean(ratio5644), dratio5644=sd(ratio5644),", 
"                              ratio4836 = mean(ratio4836), dratio4836=sd(ratio4836))", 
"        res <- rbind(res, newline)", 
"      } else {", 
"        plotwitherror(x=1/B64_56_44$L[maskB64_56_44], y=B64_56_44$DGDq2gev[maskB64_56_44], dy=B64_56_44$dDGDq2gev[maskB64_56_44],", 
"                      main=paste(\"cd sigmoid diff. decay rate , \", err, \", cd, z =\", iz, \"th\", th), ", 
"                      xlab=\"1/L\", ylab=\"24 pi^3 DGamma / Dq^2 [GeV^-3]\",", 
"                      ylim=c(min(B64_56_44$DGDq2gev[maskB64_56_44] - B64_56_44$dDGDq2gev[maskB64_56_44], ", 
"                                 B96$DGDq2gev[maskB96] - B96$dDGDq2gev[maskB96], ", 
"                                 B48$DGDq2gev[maskB48] - B48$dDGDq2gev[maskB48]),", 
"                             max(B64_56_44$DGDq2gev[maskB64_56_44] + B64_56_44$dDGDq2gev[maskB64_56_44], ", 
"                                 B96$DGDq2gev[maskB96] + B96$dDGDq2gev[maskB96],", 
"                                 B48$DGDq2gev[maskB48] + B48$dDGDq2gev[maskB48])),", 
"                      xlim=c(1/96, 1/48), xaxt=\"n\") ", 
"      }", 
"      print(c(B48$DGDq2gev[maskB48], B64_56_44$DGDq2gev[maskB64_56_44], B96$DGDq2gev[maskB96]))", 
"      axis(side=1, at=1/c(48, 64, 96), label=c(\"1/48\", \"1/64\", \"1/96\"))", 
"      legend(x=\"topleft\", legend=c(\"56, 44\", \"48, 36\"), col=c(1, 2), pch=c(1, 1), title=c(\"tsnk, tj\"))", 
"      plotwitherror(x=1/B96$L[maskB96], y=B96$DGDq2gev[maskB96], dy=B96$dDGDq2gev[maskB96],", 
"                    rep=TRUE, col=1, pch=1)", 
"      plotwitherror(x=1/B48$L[maskB48], y=B48$DGDq2gev[maskB48], dy=B48$dDGDq2gev[maskB48],", 
"                    rep=TRUE, col=2, pch=1, lwd=1.2)", 
"      plotwitherror(x=1/B64_48_36$L[maskB64_48_36], y=B64_48_36$DGDq2gev[maskB64_48_36], dy=B64_48_36$dDGDq2gev[maskB64_48_36],", 
"                    rep=TRUE, col=2, pch=1, lwd=1.2)", 
"      ", 
"    }", 
"  }", 
"}", 
sprintf("write.table(res, \"tables/DG_%s_%s_volume.csv\", row.names=FALSE)", channel, kernel), 
"```\n",
file=fileDG, append=TRUE, sep="\n")
}
}


cat("## summary slope\n", 
" We plot the slope of the fit to all four points, as well as the ratios between the differnet time separations at constant volume and different values at same time separation.\n",
file=fileDG, append=TRUE, sep="\n")

for(kernel in c("sigmoid", "erf")){
    for(channel in c("cd", "cs")) {
        cat(
        
sprintf("```{r summary_%s_%s_volume}", channel, kernel), 
"## slope",
sprintf("res <- read.table(\"tables/DG_%s_%s_volume.csv\", header=TRUE)", channel, kernel), 
"for ( err in errlist) {",
"plotwitherror(x=res$q[res$errtype==err]^2+0.005*res$iz[res$errtype==err], y=res$slope[res$errtype==err], dy=res$dslope[res$errtype==err], col=res$iz[res$errtype==err]+1, pch=res$iz[res$errtype==err]+1,", 
sprintf("              xlab=\"q^2 [GeV^2]\", ylab=\"slope\", main=paste(\"trend of infinite volume limit, %s, %s, err =\", err),", channel, kernel), 
"              ylim=range(res$slope[res$errtype==err]))", 
"lines(x=c(-1, 2), y=c(0, 0), col=\"red\")", 
"legend(title=\"iz\", col=unique(res$iz)+1, pch=unique(res$iz)+1, legend=unique(res$iz), x=\"bottomleft\")", 
"## constant volume",
"plotwitherror(x=res$q[res$errtype==err]^2+0.005*res$iz[res$errtype==err], y=res$ratio64[res$errtype==err], dy=res$dratio64[res$errtype==err], col=res$iz[res$errtype==err]+1, pch=res$iz[res$errtype==err]+1,", 
sprintf("              xlab=\"q^2 [GeV^2]\", ylab=\"ratio\", main=paste(\"ratio of B64 (56,44) / (48, 36), %s, %s, err =\", err),", channel, kernel), 
"              ylim=c(max(min(res$ratio64[res$errtype==err]), 0.5), min(max(res$ratio64[res$errtype==err]), 1.5)))", 
"lines(x=c(-1, 2), y=c(1, 1), col=\"red\")", 
"legend(title=\"iz\", col=unique(res$iz)+1, pch=unique(res$iz)+1, legend=unique(res$iz), x=\"bottomleft\")", 
"## constant time separation",
"for (z in zlist) {", 
"        plotwitherror(x=res$q[res$errtype==err & res$iz == z]^2, ",
"                      y=res$ratio5644[res$errtype==err & res$iz == z], dy=res$dratio5644[res$errtype==err & res$iz == z], ",
"                      col=1, pch=1,", 
"                      xlab=\"q^2 [GeV^2]\", ylab=\"ratio\", ",
sprintf("                     main=paste(\"ratio of (56, 44); vol 64/96 and of (48, 36); vol 48/64\n %s, %s, err =\", err, \"Z \", z),", channel, kernel), 
"                      ylim=c(max(min(res$ratio5644), 0.5), min(max(res$ratio5644), 1.5)))", 
"                      lines(x=c(-1, 2), y=c(1, 1), col=\"red\")", 
"        plotwitherror(x=res$q[res$errtype==err & res$iz == z]^2+0.005, ",
"                      y=res$ratio4836[res$errtype==err & res$iz == z], dy=res$dratio4836[res$errtype==err & res$iz == z], ",
"                      col=2, pch=2, rep=T)", 
"                      lines(x=c(-1, 2), y=c(1, 1), col=\"red\")", 
"        legend(col=c(1, 2), pch=c(1, 2), legend=c(\"(56, 44)\", \"(48, 36)\"), x=\"bottomleft\")", 
"    }",
"}",
"```\n",
file=fileDG, append=TRUE, sep="\n")
}
}


## only different time separations
fileDG <- sprintf("%s/%s-compare-stability-time-separation.Rmd", parent, 4*length(folders)+6)

cat("```{r setupsummarytimesep, include=FALSE}", 
"knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, results = FALSE)", 
"source(\"/home/gross/Documents/heavymesons/scripts_christiane/extrapolate_sigma_templates/functions_stability_plots.R\")", 
"source(\"/home/gross/Documents/heavymesons/scripts_christiane/extrapolate_sigma_templates/calc_DGDq2.R\")", 
"library(\"hadron\")", 
"zlist <- c(3, 0, 1, 2)", 
"errlist <- c(\"stat\")", 
"sumup <-FALSE", 
"```",
"", 
"# Differential Decay Rates - Volume",
file=fileDG, append=FALSE, sep="\n")


for(kernel in c("sigmoid", "erf")){
    for(channel in c("cd", "cs")) {
        for(th in c(seq(1,9), 9.5)) {
            cat(
            sprintf("## %s %s th %s", channel, kernel, th),
            sprintf("```{r timesep%s%s%s}", channel, kernel, th),
            sprintf("th1path <- \"/home/gross/Documents/heavymesons/data/%s/cB211.07.64/th%s/%s_new/outputDGammaDq2/\"", channel, th, kernel),
            sprintf("th1path2 <- \"/home/gross/Documents/heavymesons/data/%s/cB211.07.64_48_36/th%s/%s_new/outputDGammaDq2/\"", channel, th, kernel),
            "plot_stability_AA0_sets(pathlist=c(th1path, th1path2), inputfile=paste0(th1path, \"/../../DGammaDq2_sigmoid_new.in\"),", 
            "                            nnorm=3, neps=20, zlist=-1, normseq=c(1), comments=c(\"56, 44\", \"48, 36\"),",
            "                            comblist=c(0))",
            "```\n",
file=fileDG, append=TRUE, sep="\n")
        }
    }
}

## time separations and volumes
fileDG <- sprintf("%s/%s-compare-stability-volumes.Rmd", parent, 4*length(folders)+7)

cat("```{r setupsummarytimesep, include=FALSE}", 
"knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, results = FALSE)", 
"source(\"/home/gross/Documents/heavymesons/scripts_christiane/extrapolate_sigma_templates/functions_stability_plots.R\")", 
"source(\"/home/gross/Documents/heavymesons/scripts_christiane/extrapolate_sigma_templates/calc_DGDq2.R\")", 
"library(\"hadron\")", 
"zlist <- c(3, 0, 1, 2)", 
"errlist <- c(\"stat\")", 
"sumup <-FALSE", 
"```",
"", 
"# Differential Decay Rates - Volume",
file=fileDG, append=FALSE, sep="\n")


for(kernel in c("sigmoid", "erf")){
    for(channel in c("cd", "cs")) {
        for(th in c(seq(1,9), 9.5)) {
            cat(
            sprintf("## %s %s th %s", channel, kernel, th),
            sprintf("```{r timesep%s%s%s}", channel, kernel, th),
            sprintf("th1path <- \"/home/gross/Documents/heavymesons/data/%s/cB211.07.64/th%s/%s_new/outputDGammaDq2/\"", channel, th, kernel),
            sprintf("th1path2 <- \"/home/gross/Documents/heavymesons/data/%s/cB211.07.64_48_36/th%s/%s_new/outputDGammaDq2/\"", channel, th, kernel),
            sprintf("th1path3 <- \"/home/gross/Documents/heavymesons/data/%s/cB211.07.48_400/th%s/%s_new/outputDGammaDq2/\"", channel, th, kernel),
            sprintf("th1path4 <- \"/home/gross/Documents/heavymesons/data/%s/cB211.07.96/th%s/%s_new/outputDGammaDq2/\"", channel, th, kernel),
            "plot_stability_AA0_sets(pathlist=c(th1path, th1path2, th1path3, th1path4), inputfile=paste0(th1path, \"/../../DGammaDq2_sigmoid_new.in\"),", 
            "                            nnorm=3, neps=20, zlist=-1, normseq=c(1), comments=c(\"(56, 44); 64\", \"(48, 36); 64\", \"(48, 36); 48\", \"(56, 44); 96\"),",
            "                            comblist=c(0))",
            "```\n",
file=fileDG, append=TRUE, sep="\n")
        }
    }
}
        


## volume DM

#~ cat("```{r setupsummaryvolume, include=FALSE}", 
#~ "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, results = FALSE)", 
#~ "source(\"/home/gross/Documents/heavymesons/scripts_christiane/extrapolate_sigma_templates/functions_stability_plots.R\")", 
#~ "source(\"/home/gross/Documents/heavymesons/scripts_christiane/extrapolate_sigma_templates/calc_DGDq2.R\")", 
#~ "library(\"hadron\")", 
#~ "zlist <- c(4, 0, 1, 2, 3)", 
#~ "errlist <- c(\"stat\")", 
#~ "sumup <-FALSE", 
#~ "```",
#~ "", 
#~ "# Differential Decay Rates - Volume",
#~ file=fileDM, append=FALSE, sep="\n")


#~ for(kernel in c("sigmoid", "erf")){
#~     for(channel in c("cd", "cs")) {
#~         cat(sprintf("## %s %s\n", channel, kernel), 
#~ sprintf("```{r decayrate_%s_volume_%s}", channel, kernel), 
#~ sprintf("B64 <- read.table(\"tables/DM_complete_cB64_converted_%s_finvol_%s.csv\")", channel, kernel), 
#~ sprintf("B96 <- read.table(\"tables/DM_complete_cB96_converted_%s_finvol_%s.csv\")", channel, kernel), 
#~ sprintf("B48 <- read.table(\"tables/DM_complete_cB48_converted_%s_finvol_%s.csv\")", channel, kernel), 
#~ "", 
#~ "", 
#~ sprintf("B64dat <- readRDS(\"tables/DM_complete_cB64_converted_%s_finvol_%s.RDS\")", channel, kernel), 
#~ sprintf("B96dat <- readRDS(\"tables/DM_complete_cB96_converted_%s_finvol_%s.RDS\")", channel, kernel), 
#~ sprintf("B48dat <- readRDS(\"tables/DM_complete_cB48_converted_%s_finvol_%s.RDS\")", channel, kernel), 
#~ "", 
#~ "B64 <- na.omit(B64)", 
#~ "B96 <- na.omit(B96)", 
#~ "B48 <- na.omit(B48)", 
#~ "", 
#~ "B64 <- B64[order(B64$iz, B64$errtype, B64$icomb, B64$q), ]", 
#~ "B96 <- B96[order(B96$iz, B96$errtype, B96$icomb, B96$q), ]", 
#~ "B48 <- B48[order(B48$iz, B48$errtype, B48$icomb, B48$q), ]", 
#~ "", 
#~ "B64$L <- rep(64, length(B64$iz))", 
#~ "B96$L <- rep(96, length(B96$iz))", 
#~ "B48$L <- rep(48, length(B48$iz))", 
#~ "", 
#~ "", 
#~ "tablelist <- list(B48, B64, B96)", 
#~ "enslist <- list(B48dat, B64dat, B96dat)", 
#~ "", 
#~ "", 
#~ "", 
#~ "fnlin <- function(par, x, boot.R, ...) par[1] + par[2] * x", 
#~ "fncon <- function(par, x, boot.R, ...) par[1]", 
#~ "", 
#~ "fitfn <- fnlin", 
#~ "par.guess <- c(1, 1)", 
#~ "comment <- \"linear\"", 
#~ "", 
#~ "res <- data.frame(th=c(), errtype=c(), iz=c(), lim=c(), dlim=c(), sloipe=c(), dslope=c(), chi=c(), q=c())", 
#~ "", 
#~ "", 
#~ "pcol <- col2rgb(\"gray\", alpha=TRUE)/255 ", 
#~ "pcol[4] <- 0.65", 
#~ "pcol <- rgb(red=pcol[1],green=pcol[2],blue=pcol[3],alpha=pcol[4])", 
#~ "", 
#~ "confcounter <- 1", 
#~ "", 
#~ "", 
#~ "", 
#~ "for(th in c(seq(1, 9), 9.5)) {", 
#~ "  for(iz in zlist) {", 
#~ "    for (err in errlist){", 
#~ "      ", 
#~ "      title <- sprintf(\"th%s q^2 %.3f GeV^2 iz %d errtype %s\", th, B64$q[abs(B64$th - th) < 1e-2][1]^2, iz, err)", 
#~ "      print(title)", 
#~ "      bsamples <- array(NA, dim=c(1000, length(enslist)))", 
#~ "      y <- c()", 
#~ "      dy <- c()", 
#~ "      x <- c()", 
#~ "      for (index in seq(1, length(enslist))) {", 
#~ "        if(!is.na(tablelist[[index]]$q[abs(tablelist[[index]]$th - th) < 1e-2][1]^2)) {", 
#~ "          mask <- abs(enslist[[index]]$th - th) < 1e-2 & enslist[[index]]$iz == iz & enslist[[index]]$errtype==err", 
#~ "          masktable <- abs(tablelist[[index]]$th - th) < 1e-2 & tablelist[[index]]$iz == iz & tablelist[[index]]$errtype==err", 
#~ "          mask[is.na(mask)] <- F", 
#~ "          y[index] <- enslist[[index]]$DGDq2gev[mask]", 
#~ "          dy[index] <- enslist[[index]]$dDGDq2gev[mask]", 
#~ "          try(bsamples[, index] <- enslist[[index]]$datgev[, mask])", 
#~ "          x[index] <- tablelist[[index]]$L[1]", 
#~ "        }", 
#~ "      }", 
#~ "      ", 
#~ "      maskB64 <- B64$iz==iz & B64$icomb==0 & B64$errtype==err & B64$th==th", 
#~ "      maskB96 <- B96$iz==iz & B96$icomb==0 & B96$errtype==err & B96$th==th", 
#~ "      maskB48 <- B48$iz==iz & B48$icomb==0 & B48$errtype==err & B48$th==th", 
#~ "      ", 
#~ "      fitresult <- try(bootstrap.nlsfit(x=1/x, y=y, bsamples=bsamples, fn=fitfn, par.guess=par.guess))", 
#~ "      ", 
#~ "      if(!inherits(x = fitresult, what=\"try-error\")){", 
#~ "        plot(fitresult,", 
#~ sprintf("             main=paste(\"%s %s diff. decay rate , \", err, \", z =\", iz, \"th\", th), ", channel, kernel), 
#~ "             xlab=\"1/L\", ylab=\"24 pi^3 DGamma / Dq^2 [GeV^-3]\", xaxt=\"n\")", 
#~ "        plotwitherror(x=1/B64$L[maskB64], y=B64$DGDq2gev[maskB64], dy=B64$dDGDq2gev[maskB64], rep=T, pch=1)", 
#~ "        newline <- data.frame(th=th, errtype=err, iz=iz, ", 
#~ "                              lim=fitresult$t0[1], dlim=fitresult$se[1], ", 
#~ "                              slope=fitresult$t0[2], dslope=fitresult$se[2], ", 
#~ "                              chi=fitresult$chisqr / fitresult$dof, q=B64$q[abs(B64$th - th) < 1e-2][1])", 
#~ "        res <- rbind(res, newline)", 
#~ "      } else {", 
#~ "        plotwitherror(x=1/B64$L[maskB64], y=B64$DGDq2gev[maskB64], dy=B64$dDGDq2gev[maskB64],", 
#~ "                      main=paste(\"cd sigmoid diff. decay rate , \", err, \", cd, z =\", iz, \"th\", th), ", 
#~ "                      xlab=\"1/L\", ylab=\"24 pi^3 DGamma / Dq^2 [GeV^-3]\",", 
#~ "                      ylim=c(min(B64$DGDq2gev[maskB64] - B64$dDGDq2gev[maskB64], ", 
#~ "                                 B96$DGDq2gev[maskB96] - B96$dDGDq2gev[maskB96], ", 
#~ "                                 B48$DGDq2gev[maskB48] - B48$dDGDq2gev[maskB48]),", 
#~ "                             max(B64$DGDq2gev[maskB64] + B64$dDGDq2gev[maskB64], ", 
#~ "                                 B96$DGDq2gev[maskB96] + B96$dDGDq2gev[maskB96],", 
#~ "                                 B48$DGDq2gev[maskB48] + B48$dDGDq2gev[maskB48])),", 
#~ "                      xlim=c(1/96, 1/48), xaxt=\"n\") ", 
#~ "      }", 
#~ "      print(c(B48$DGDq2gev[maskB48], B64$DGDq2gev[maskB64], B96$DGDq2gev[maskB96]))", 
#~ "      axis(side=1, at=1/c(48, 64, 96), label=c(\"1/48\", \"1/64\", \"1/96\"))", 
#~ "      legend(x=\"topleft\", legend=c(\"56, 44\", \"48, 32\"), col=c(1, 2), pch=c(1, 1), title=c(\"tsnk, tj\"))", 
#~ "      plotwitherror(x=1/B96$L[maskB96], y=B96$DGDq2gev[maskB96], dy=B96$dDGDq2gev[maskB96],", 
#~ "                    rep=TRUE, col=1, pch=1)", 
#~ "      plotwitherror(x=1/B48$L[maskB48], y=B48$DGDq2gev[maskB48], dy=B48$dDGDq2gev[maskB48],", 
#~ "                    rep=TRUE, col=2, pch=1)", 
#~ "      ", 
#~ "    }", 
#~ "  }", 
#~ "}", 
#~ sprintf("write.table(res, \"tables/DG_%s_%s_volume.csv\", row.names=FALSE)", channel, kernel), 
#~ "```\n",
#~ file=fileDM, append=TRUE, sep="\n")
#~ }
#~ }


#~ cat("## summary slope\n", file=fileDM, append=TRUE, sep="\n")

#~ for(kernel in c("sigmoid", "erf")){
#~     for(channel in c("cd", "cs")) {
#~         cat(
#~ sprintf("```{r summary_%s_%s_volume}", channel, kernel), 
#~ sprintf("res <- read.table(\"tables/DG_%s_%s_volume.csv\", header=TRUE)", channel, kernel), 
#~ "plotwitherror(x=res$q^2+0.005*res$iz, y=res$slope, dy=res$dslope, col=res$iz+1, pch=res$iz+1,", 
#~ sprintf("              xlab=\"q^2 [GeV^2]\", ylab=\"slope\", main=\"trend of infinite volume limit, %s, %s\",", channel, kernel), 
#~ "              ylim=range(res$slope))", 
#~ "lines(x=c(-1, 2), y=c(0, 0), col=\"red\")", 
#~ "legend(title=\"iz\", col=unique(res$iz)+1, pch=unique(res$iz)+1, legend=unique(res$iz), x=\"bottomleft\")", 
#~ "```\n",
#~ file=fileDM, append=TRUE, sep="\n")
#~ }
#~ }

#~ ## TODO: write continuum limit
