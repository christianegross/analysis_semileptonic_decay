parent <- "~/Documents/heavymesons/data/newinput"
folders <- c("cB211.07.64", "cB211.07.96", "cC211.06.80", "cD211.054.96", "cE211.044.112", "cB211.07.48_300", "cB211.07.48_400", "cB211.07.64_48_36")
nameshort <- c("B64", "B96", "C80", "D96", "E112", "B48_300", "B48_400", "B64_48_36")
subfolders  <- c("th1", "th2", "th3", "th4", "th5", "th6", "th7", "th8", "th9", "th9.5")

infos <- read.table("parameters_input_files.csv", row.names=1, header=TRUE, sep=",", comment.char = "#")
#~ print(infos)
names(infos)



## continuum limit DG


fileDG <- sprintf("%s/%s-integral-DGDq2.Rmd", parent, 4*length(folders)+8)
fileDM <- sprintf("%s/%s-integral-DMDq2.Rmd", parent, 4*length(folders)+8)

## overview DG
cat("```{r setupsummaryintegral, include=FALSE}", 
"knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, results = FALSE)", 
"source(\"/home/gross/Documents/heavymesons/scripts_christiane/extrapolate_sigma_templates/functions_stability_plots.R\")", 
"source(\"/home/gross/Documents/heavymesons/scripts_christiane/extrapolate_sigma_templates/calc_DGDq2.R\")", 
"source(\"/home/gross/Documents/heavymesons/scripts_christiane/extrapolate_sigma_templates/integrate.R\")", 
"source(\"/home/gross/Documents/heavymesons/scripts_christiane/integrals/splineintegration_functions.R\")", 
"library(\"hadron\")", 
"zlist <- c(3, 0, 1, 2)", 
"errlist <- c(\"stat\", \"sys\", \"vol\", \"tot\")", 
"doplot <- T", 
"fnlin <- function(par, x, boot.R, ...) par[1] + par[2] * x", 
"fncon <- function(par, x, boot.R, ...) par[1]", 
"", 
"```",
"", 
"# Integral",
file=fileDG, append=FALSE, sep="\n")


for(kernel in c("sigmoid", "erf", "aic")){
    for(channel in c("cd", "cs")) {
        cat(sprintf("## %s %s\n", channel, kernel),
sprintf("```{r integral%s%s}", channel, kernel),
sprintf("B64 <- readRDS(\"tables/DG_complete_cB64_converted_%s_finvol_%s.RDS\")", channel, kernel), 
sprintf("C80 <- readRDS(\"tables/DG_complete_cC80_converted_%s_finvol_%s.RDS\")", channel, kernel), 
sprintf("D96 <- readRDS(\"tables/DG_complete_cD96_converted_%s_finvol_%s.RDS\")", channel, kernel), 
sprintf("E112 <- readRDS(\"tables/DG_complete_cE112_converted_%s_finvol_%s.RDS\")", channel, kernel), 
sprintf("contlim <- readRDS(\"tables/DG_contlimit_%s_%s_linear.RDS\")", channel, kernel),
"", 
sprintf("B64table <- read.table(\"tables/DG_complete_cB64_converted_%s_finvol_%s.csv\", header=TRUE)", channel, kernel), 
sprintf("C80table <- read.table(\"tables/DG_complete_cC80_converted_%s_finvol_%s.csv\", header=TRUE)", channel, kernel), 
sprintf("D96table <- read.table(\"tables/DG_complete_cD96_converted_%s_finvol_%s.csv\", header=TRUE)", channel, kernel), 
sprintf("E112table <- read.table(\"tables/DG_complete_cE112_converted_%s_finvol_%s.csv\", header=TRUE)", channel, kernel),  
sprintf("contlimtable <- read.table(\"tables/DG_contlimit_%s_%s_linear.csv\")", channel, kernel),
"", 
"B64$L <- rep(64, length(B64$Nt))", 
"C80$L <- rep(80, length(C80$Nt))", 
"D96$L <- rep(96, length(D96$Nt))", 
"E112$L <- rep(112, length(E112$Nt))", 
"", 
"", 
"enslist <- list(B64, C80, D96, E112, contlim)", 
"tablelist <- list(B64table, C80table, D96table, E112table, contlimtable)", 
"names <- c(\"B64\", \"C80\", \"D96\", \"E112\", \"contlim\")", 
"", 
"", 
"res <- data.frame(name=NA, ensno=NA, errtype=NA, iz=NA, intspline=NA, dintspline=NA, intbsspline=NA, inttrap=NA, dinttrap=NA, intbstrap=NA)", 
"reslist <- list()", 
"upperboundarycd <- 0.8724",
"upperboundarycs <- 0.7770",
"", 
"",
"confcounter <- 1", 
"thetas <- c(seq(1, 9), 9.5)", 
"for (i in seq_along(enslist)){", 
"  ens <- enslist[[i]]", 
"  mytable <- tablelist[[i]]", 
"  for (errtype in errlist) {", 
"    for (iz in zlist) {", 
"      title <- sprintf(\"ens %s iz %d errtype %s\", names[i], iz, errtype)", 
"      titlelist <- sprintf(\"ens%s-iz%d-errtype%s\", names[i], iz, errtype)", 
"      print(title)", 
"      ## we include the error on a in the continuum limit", 
"      bsamples <- array(rep(0, 11000), dim=c(1000, 11))", 
"      y <- c(0)", 
"      dy <- c(0)", 
"      x <- c(0)", 
"      dx <- c(0)", 
"      for (index in seq_along(thetas)) {", 
"        th <- thetas[index]", 
"        if(!is.na(mytable$q[abs(mytable$th - th) < 1e-2][1]^2) && names[i] != \"contlim\") {", 
"          mask <- abs(ens$th - th) < 1e-2 & ens$iz == iz & ens$errtype==errtype", 
"          mask[is.na(mask)] <- F", 
"          masktable <- abs(mytable$th - th) < 1e-2 & mytable$iz == iz & mytable$errtype==errtype", 
"          masktable[is.na(masktable)] <- F", 
"          y[index+1] <- ens$DGDq2gev[mask]", 
"          dy[index+1] <- ens$dDGDq2gev[mask]", 
"          try(bsamples[, index+1] <- ens$datgev[, mask])", 
"          x[index+1] <- mytable$q[masktable]^2", 
"        } else if(!is.na(mytable$q[abs(mytable$th - th) < 1e-2][1]^2) && names[i] == \"contlim\") {
          mask <- abs(ens$th - th) < 1e-2 & ens$iz == iz & ens$errtype==errtype
          mask[is.na(mask)] <- F
          masktable <- abs(mytable$th - th) < 1e-2 & mytable$iz == iz & mytable$errtype==errtype
          masktable[is.na(masktable)] <- F
          y[index+1] <- ens$fit[[which(mask)]]$t0
          dy[index+1] <- ens$fit[[which(mask)]]$se
          try(bsamples[, index+1] <- ens$fit[[which(mask)]]$t[, 1])
          x[index+1] <- mytable$q[masktable]^2
        } else {", 
"          print(sprintf(\"failure for th %s\", index))", 
"        }", 
"      }", 
"      ",
file=fileDG, append=TRUE, sep="\n") 
if(channel=="cs") {
    cat(
"      spline <- interpSpline(x, y)", 
"      len <- length(x)", 
"      slopebs <- (bsamples[, len] - bsamples[, len-1])/(x[len] - x[len-1])", 
"      yupperbs <- bsamples[, len] + (upperboundarycs-x[len]) * slopebs", 
"      ", 
"      plotwitherror(x=x, y=y, dy=dy, main=title, xlab=\"q^2\", ylab=\"DGammaDq^2\")", 
"      xval <- seq(min(x), upperboundarycs, length.out=500)", 
"      lines(x=xval, y=predict(object=spline, x=xval)$y, col=\"red\", lty=2)", 
"      lines(x[1:(len-1)], y[1:(len-1)], col=\"blue\", lty=3)", 
"      # plotwitherror(x=upperboundarycs, y=mean(yupperbs), dy=sd(yupperbs), col=\"blue\", rep=T)", 
"      lines(x=c(x[len-1], upperboundarycs), y=c(y[len-1], mean(yupperbs)), col=\"blue\", lty=3)", 
"      plotwitherror(x=x, y=y, dy=dy, rep=TRUE)", 
"      legend(x=\"topright\", legend=c(\"meas\", \"spline\", \"trapezoidal\"), col=c(1, \"red\", \"blue\"), pch=c(1, NA, NA), lty=c(NA, 2, 3))", 
"      meanintspline <- splineintegral(yval=y, xval=x, replacelower = TRUE, replaceindex = 10, lowerlimit = upperboundarycs)", 
"      meaninttrap <- trapezoidal(yval=y, xval=x, replacelower = TRUE, replaceindex = 10, lowerlimit = upperboundarycs)", 
"      bsintspline <- apply(X=bsamples, MARGIN=1, FUN=splineintegral, xval=x, replacelower = TRUE, replaceindex = 10, lowerlimit = upperboundarycs)", 
"      bsinttrap <- apply(X=bsamples, MARGIN=1, FUN=trapezoidal, xval=x, replacelower = TRUE, replaceindex = 10, lowerlimit = upperboundarycs)", 
"      res <- rbind(res, data.frame(name=names[i], ensno=i, errtype=errtype, iz=iz, ", 
"                                   intspline=meanintspline, dintspline=sd(bsintspline), intbsspline=mean(bsintspline), ", 
"                                   inttrap=meaninttrap, dinttrap=sd(bsinttrap), intbstrap=mean(bsinttrap)))",
file=fileDG, append=TRUE, sep="\n")
}
if(channel=="cd") {
    cat(
"      spline <- interpSpline(x, y)", 
"      len <- length(x)", 
"      slopebs <- (bsamples[, len] - bsamples[, len-1])/(x[len] - x[len-1])", 
"      yupperbs <- bsamples[, len] + (upperboundarycd-x[len]) * slopebs", 
"      ", 
"      plotwitherror(x=x, y=y, dy=dy, main=title, xlab=\"q^2\", ylab=\"DGammaDq^2\")", 
"      xval <- seq(min(x), upperboundarycd, length.out=500)", 
"      lines(x=xval, y=predict(object=spline, x=xval)$y, col=\"red\", lty=2)", 
"      lines(x, y, col=\"blue\", lty=3)", 
"      # plotwitherror(x=upperboundarycd, y=mean(yupperbs), dy=sd(yupperbs), col=\"blue\", rep=T)", 
"      lines(x=c(x[len], upperboundarycd), y=c(y[len], mean(yupperbs)), col=\"blue\", lty=3)", 
"      plotwitherror(x=x, y=y, dy=dy, rep=TRUE)", 
"      legend(x=\"topright\", legend=c(\"meas\", \"spline\", \"trapezoidal\"), col=c(1, \"red\", \"blue\"), pch=c(1, NA, NA), lty=c(NA, 2, 3))", 
"      ", 
"      meanintspline <- splineintegral(yval=y, xval=x, continue = T, higherlimit = upperboundarycd)", 
"      meaninttrap <- trapezoidal(yval=y, xval=x, continue = T, higherlimit = upperboundarycd)", 
"      bsintspline <- apply(X=bsamples, MARGIN=1, FUN=splineintegral, xval=x, continue = T, higherlimit = upperboundarycd)", 
"      bsinttrap <- apply(X=bsamples, MARGIN=1, FUN=trapezoidal, xval=x, continue = T, higherlimit = upperboundarycd)", 
"      res <- rbind(res, data.frame(name=names[i], ensno=i, errtype=errtype, iz=iz, ", 
"                                   intspline=meanintspline, dintspline=sd(bsintspline), intbsspline=mean(bsintspline), ", 
"                                   inttrap=meaninttrap, dinttrap=sd(bsinttrap), intbstrap=mean(bsinttrap)))",
file=fileDG, append=TRUE, sep="\n")
}
cat(
"      ",
"      reslist[[paste0(names[i], errtype, \"iz\", iz, \"bsintspline\")]] <- bsintspline",
"      reslist[[paste0(names[i], errtype, \"iz\", iz, \"spline\")]] <- spline",
"      reslist[[paste0(names[i], errtype, \"iz\", iz, \"bsinttrap\")]] <- bsinttrap",
"      confcounter <- confcounter + 1", 
"    }", 
"  }", 
"}", 
"", 
"res <- res[-1, ]", 
"res",
sprintf("write.table(x=res, file=\"tables/DG_%s_%s_int.csv\", row.names=F, col.names=T)", channel, kernel),
"",
"reslist$info <- res",
sprintf("saveRDS(object=reslist, file=\"tables/DG_%s_%s_int.RDS\") ", channel, kernel),
"", 
"```\n",
"\n\n### Summary\n\n",
"```{r}\n",
"for(errtype in errlist) {", 
"  mask <- res$errtype==errtype", 
"  plotwitherror(x=res$iz[mask]+ 0.1*(res$ensno[mask]-1), y=res$intspline[mask], dy=res$dintspline[mask], col=res$ensno[mask], pch=res$ensno[mask], ", 
"                main=paste(\"err =\", errtype), xlab=\"Z\", ylab=\"Gamma [a.u.]\")", 
"  legend(legend=names, x=\"topleft\", col=seq(1, max(res$ensno)), pch=seq(1, max(res$ensno)))",
"}", 
"```\n",
file=fileDG, append=TRUE, sep="\n")
}
}

folders <- c("cB211.07.64")
nameshort <- c("B64")
subfolders  <- c("th2_su", "th4_su", "th6_su", "th8_su", "th9.5_su", "th1")


for(kernel in c("sigmoid", "erf")){
    for(channel in c("su")) {
        cat(sprintf("## %s %s\n", channel, kernel),
sprintf("```{r integral%s%s}", channel, kernel),
"errlist <- c(\"stat\", \"sys\")",
sprintf("B64 <- readRDS(\"tables/DG_complete_cB64_converted_%s_finvol_%s.RDS\")", channel, kernel), 
"", 
sprintf("B64table <- read.table(\"tables/DG_complete_cB64_converted_%s_finvol_%s.csv\", header=TRUE)", channel, kernel), 
"", 
"B64$L <- rep(64, length(B64$Nt))", 
"", 
"", 
"enslist <- list(B64)", 
"tablelist <- list(B64table)", 
"names <- c(\"B64\")", 
"upperboundarysu <- 0.0069",
"", 
"", 
"res <- data.frame(name=NA, ensno=NA, errtype=NA, iz=NA, intspline=NA, dintspline=NA, intbsspline=NA, inttrap=NA, dinttrap=NA, intbstrap=NA)", 
"reslist <- list()", 
"", 
"", 
"", 
"confcounter <- 1", 
"thetas <- c(\"2_su\", \"4_su\", \"6_su\", \"8_su\", \"9.5_su\", \"1\")", 
"ntheta <- length(thetas)+1",
"for (i in seq_along(enslist)){", 
"  ens <- enslist[[i]]", 
"  mytable <- tablelist[[i]]", 
"  for (errtype in errlist) {", 
"    for (iz in zlist) {", 
"      title <- sprintf(\"ens %s iz %d errtype %s\", names[i], iz, errtype)", 
"      titlelist <- sprintf(\"ens%s-iz%d-errtype%s\", names[i], iz, errtype)", 
"      print(title)", 
"      ## we include the error on a in the continuum limit", 
"      bsamples <- array(rep(0, ntheta*1000), dim=c(1000, ntheta))", 
"      y <- c(0)", 
"      dy <- c(0)", 
"      x <- c(0)", 
"      dx <- c(0)", 
"      for (index in seq_along(thetas)) {", 
"        th <- thetas[index]", 
"        if(!is.na(mytable$q[as.character(mytable$th) == as.character(th)][1]^2) && names[i] != \"contlim\") {", 
"          mask <- as.character(ens$th) == as.character(th) & ens$iz == iz & ens$errtype==errtype", 
"          mask[is.na(mask)] <- F", 
"          masktable <- as.character(mytable$th) == as.character(th) & mytable$iz == iz & mytable$errtype==errtype", 
"          masktable[is.na(masktable)] <- F", 
"          y[index+1] <- ens$DGDq2gev[mask]", 
"          dy[index+1] <- ens$dDGDq2gev[mask]", 
"          try(bsamples[, index+1] <- ens$datgev[, mask])", 
"          x[index+1] <- mytable$q[masktable]^2", 
"        } else if(!is.na(mytable$q[as.character(mytable$th) == as.character(th)][1]^2) && names[i] == \"contlim\") {
          mask <- as.character(ens$th) == as.character(th) & ens$iz == iz & ens$errtype==errtype
          mask[is.na(mask)] <- F
          masktable <- as.character(mytable$th) == as.character(th) & mytable$iz == iz & mytable$errtype==errtype
          masktable[is.na(masktable)] <- F
          y[index+1] <- ens$fit[[which(mask)]]$t0
          dy[index+1] <- ens$fit[[which(mask)]]$se
          try(bsamples[, index+1] <- ens$fit[[which(mask)]]$t[, 1])
          x[index+1] <- mytable$q[masktable]^2
        } else {", 
"          print(sprintf(\"failure for th %s\", index))", 
"        }", 
"      }", 
"      ", 
"      spline <- interpSpline(x, y)", 
"      len <- length(x)", 
"      slopebs <- (bsamples[, len] - bsamples[, len-1])/(x[len] - x[len-1])", 
"      yupperbs <- bsamples[, len] + (upperboundarysu-x[len]) * slopebs", 
"      ", 
"      plotwitherror(x=x, y=y, dy=dy, main=title, xlab=\"q^2\", ylab=\"DGammaDq^2\")", 
"      xval <- seq(min(x), upperboundarysu, length.out=500)", 
"      lines(x=xval, y=predict(object=spline, x=xval)$y, col=\"red\", lty=2)", 
"      lines(x[1:(len-1)], y[1:(len-1)], col=\"blue\", lty=3)", 
"      lines(x=c(x[len-1], upperboundarysu), y=c(y[len-1], mean(yupperbs)), col=\"blue\", lty=3)", 
"      plotwitherror(x=x, y=y, dy=dy, rep=TRUE)", 
"      legend(x=\"topright\", legend=c(\"meas\", \"spline\", \"trapezoidal\"), col=c(1, \"red\", \"blue\"), pch=c(1, NA, NA), lty=c(NA, 2, 3))", 
"      meanintspline <- splineintegral(yval=y, xval=x, replacelower = TRUE, replaceindex = 6, lowerlimit = upperboundarysu)", 
"      meaninttrap <- trapezoidal(yval=y, xval=x, replacelower = TRUE, replaceindex = 6, lowerlimit = upperboundarysu)", 
"      bsintspline <- apply(X=bsamples, MARGIN=1, FUN=splineintegral, xval=x, replacelower = TRUE, replaceindex = 6, lowerlimit = upperboundarysu)", 
"      bsinttrap <- apply(X=bsamples, MARGIN=1, FUN=trapezoidal, xval=x, replacelower = TRUE, replaceindex = 6, lowerlimit = upperboundarysu)", 
"      res <- rbind(res, data.frame(name=names[i], ensno=i, errtype=errtype, iz=iz, ", 
"                                   intspline=meanintspline, dintspline=sd(bsintspline), intbsspline=mean(bsintspline), ", 
"                                   inttrap=meaninttrap, dinttrap=sd(bsinttrap), intbstrap=mean(bsinttrap)))", 
"      reslist[[paste0(names[i], errtype, \"iz\", iz, \"bsintspline\")]] <- bsintspline",
"      reslist[[paste0(names[i], errtype, \"iz\", iz, \"spline\")]] <- spline",
"      reslist[[paste0(names[i], errtype, \"iz\", iz, \"bsinttrap\")]] <- bsinttrap",
"      ",
"      ", 
"      confcounter <- confcounter + 1", 
"    }", 
"  }", 
"}", 
"", 
"res <- res[-1, ]", 
"res",
sprintf("write.table(x=res, file=\"tables/DG_%s_%s_int.csv\", row.names=F, col.names=T)", channel, kernel),
"",
"reslist$info <- res",
sprintf("saveRDS(object=reslist, file=\"tables/DG_%s_%s_int.RDS\") ", channel, kernel),
"", 
"```\n",
"\n\n### Summary\n\n",
"```{r}\n",
"for(errtype in errlist) {", 
"  mask <- res$errtype==errtype", 
"  plotwitherror(x=res$iz[mask]+ 0.1*(res$ensno[mask]-1), y=res$intspline[mask], dy=res$dintspline[mask], col=res$ensno[mask], pch=res$ensno[mask], ", 
"                main=paste(\"err =\", errtype), xlab=\"Z\", ylab=\"Gamma [a.u.]\")", 
"  legend(legend=names, x=\"topleft\", col=seq(1, max(res$ensno)), pch=seq(1, max(res$ensno)))",
"}", 
"```\n",
file=fileDG, append=TRUE, sep="\n")
}
}
