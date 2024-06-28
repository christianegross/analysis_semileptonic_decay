library("hadron")

## cdf needed for determining mean, intervals with AIC
cdf <- function(y, means, sds, weights) {
  num <- sum(weights * pnorm(y, means, sds))
  den <- sum(weights)
  return(num / den)
}

## cdf - quantile, needed for uniroot

AICquantile <- function(y, means, sds, weights, quantile) cdf(y, means, sds, weights) - quantile

## apply requires the argument to be iterated over in first place
getquantile <- function(means, sds, weights, quantile, interval) {
    root <- try(uniroot(f = AICquantile, quantile = quantile, means = means, sds = sds,
        weights = weights, interval = interval, tol = 1e-12))
    if (inherits(root, "try-error")) {
        res <- NA
    } else {
        res <- root$root
    }
    return(res)
}


## do fits with different functions and determine weights with akaike information criterion
## The parameter we want to determine is the first parameter of all fitfunctions, in this case the constant of a polynomial extrapolating to zero
combinefunctionsaic <- function(xlist, ylist, bsampleslist, functions, masklist, pars, verbose = TRUE) {
    stopifnot (length(functions) == length(pars))
    stopifnot (length(xlist) == length(ylist))
    boot.R <- length(bsampleslist[[1]][, 1])
    
    ncombs <- length(functions)
    weights <- rep(NA, ncombs)
    means <- rep(NA, ncombs)
    sds <- rep(NA, ncombs)
    bootstraps <- array(NA, dim = c(boot.R, ncombs))
    savefits <- data.frame(index = c(),
                            mean = c(), sd = c(), chisqr = c(), p = c(), weight = c())
    fitresults <- list()
                            
    for (index in seq(1, length(functions))) {
            x <- xlist[[index]]
            y <- ylist[[index]]
            bsamples <- bsampleslist[[index]]
            stopifnot (length(x) == length(y))
#~             plotwitherror(x=x, y=y, dy=apply(X=bsamples, FUN=sd, MARGIN=2))
#~             print(functions[[index]])
#~             print(pars[[index]])
#~             print(which(masklist[[index]]))
            tmp <- try(bootstrap.nlsfit(x=x, y=y, bsamples=bsamples, 
                        fn=functions[[index]], par.guess=pars[[index]],
                        mask=masklist[[index]]))#,
                        #CovMatrix=NULL))
#~             print(summary(tmp))
            if (!inherits(tmp, "try-error")) {
                weights[index] <- exp(-0.5 * (tmp$chisqr + length(pars[[index]]) - length(x)))
                means[index] <- tmp$t0[1]
                sds[index] <- tmp$se[1]
                bootstraps[, index] <- tmp$t[, 1]
                newline <- data.frame(index = index, 
                            mean = tmp$t0[1], sd = tmp$se[1],
                            chisqr = tmp$chisqr, p = tmp$Qval, weight = weights[index])
                savefits <- rbind(savefits, newline)
                fitresults[[index]] <- tmp
            } else {
                weights[index] <- 0
                means[index] <- 0
                sds[index] <- 0
                bootstraps[, index] <- rep(0, length(bsamples[, 1]))
                cat("there was a problem with index ", index, "\n")
            }
    }

    # failsafe if no fit succeeds
    if(length(savefits$index) == 0) return(list(FALSE, FALSE, FALSE, FALSE))

    savefits$weight <- savefits$weight / sum(savefits$weight)
    weights <- weights / sum(weights)
    
    
    ## determine mean and error as median and 16%/84% quantile of the cdf
    interval <- c(0.9 * min(bootstraps), 1.1 * max(bootstraps))

    cdfresmean <- uniroot(f = AICquantile, quantile = 0.5, means = means, sds = sds,
        weights = weights, interval = interval, tol = 1e-12)
    cdfres16 <- uniroot(f = AICquantile, quantile = 0.16, means = means, sds = sds,
        weights = weights, interval = interval, tol = 1e-12)
    cdfres84 <- uniroot(f = AICquantile, quantile = 0.84, means = means, sds = sds,
        weights = weights, interval = interval, tol = 1e-12)

    cdfresmeanerr <- abs(0.5 * (cdfres84$root - cdfres16$root))

    cdfbootmeanbs <- apply(X = bootstraps, MARGIN = 1, FUN = getquantile, sds = sds,
        weights = weights, quantile = 0.5, interval = interval)
    cdfboot16 <- apply(X = bootstraps, MARGIN = 1, FUN = getquantile, sds = sds,
        weights = weights, quantile = 0.16, interval = interval)
    cdfboot84 <- apply(X = bootstraps, MARGIN = 1, FUN = getquantile, sds = sds,
        weights = weights, quantile = 0.84, interval = interval)

    cdfbootmean <- mean(cdfbootmeanbs, na.rm = T)
    cdfbooterrstat <- sd(cdfbootmeanbs, na.rm = T)
    cdfbooterrtot <- abs(0.5 * (mean(cdfboot84, na.rm = T) - mean(cdfboot16, na.rm = T)))

    if (verbose) {
        print(savefits)
        print(paste("min = ", cdf(interval[1], means = means, sds = sds, weights = weights),
        " max = ", cdf(interval[2], means = means, sds = sds, weights = weights)))
        print(data.frame(t0 = cdfresmean$root, setot = cdfresmeanerr, 
                mean = cdfbootmean, errstat = cdfbooterrstat,
                booterrtot = cdfbooterrtot, booterrsys = sqrt(cdfbooterrtot^2 - cdfbooterrstat^2),
                bias = cdfresmean$root - cdfbootmean, m16 = mean(cdfboot16), m84 = mean(cdfboot84)))
    }
    
    ## determine median and error by weighting the results
    
    wmean <- sum(weights*means)
    werr <- sqrt(sum(weights*sds^2) + sum(weights*(means-wmean)^2))
    wbootstraps <- apply(MARGIN=1, FUN=sum, X=bootstraps*rep(weights, each=length(bootstraps[, 1])))

    return(list(fitresults=fitresults, 
                savefits=savefits,
                cdfsummary=data.frame(t0 = cdfresmean$root, setot = cdfresmeanerr, 
                    mean = cdfbootmean, errstat = cdfbooterrstat,
                    booterrtot = cdfbooterrtot, booterrsys = sqrt(cdfbooterrtot^2 - cdfbooterrstat^2),
                    bias = cdfresmean$root - cdfbootmean, m16 = mean(cdfboot16), m84 = mean(cdfboot84)),
                cdfbootsamples = list(bootmean=cdfbootmeanbs, boot16=cdfboot16, boot84=cdfboot84),
                weights=weights,
                weighted=list(mean=wmean, sd=werr, bootsamples=wbootstraps)
        ))


}

plotlistfits <- function(fitlist, xlim=c(NA), ylim=c(NA), cols=seq(1, length(fitlist)), 
                         drawlegend=F, legendargs=list(), legendwide=F, ...) {
    if (is.na(xlim[1])) xlim <- minmax(unlist(sapply(X=fitlist, FUN=getElement, name="x"), use.names=F))
    if (is.na(ylim[1])) ylim <- minmax(unlist(sapply(X=fitlist, FUN=getElement, name="y"), use.names=F))
    if(drawlegend) {
        xwidth <- xlim[2] - xlim[1]
        ywidth <- ylim[2] - ylim[1]
        if(legendwide) ylim[2] <- ylim[2] + 0.1*ywidth
        if(!legendwide) xlim[2] <- xlim[2] + 0.1*xwidth
    }
    plot(NA, xlim=xlim, ylim=ylim, ...)
    xseq <- seq(from=xlim[1]*1-0.1*sign(xlim[1]), to=xlim[2]*1+0.1*sign(xlim[2]), length.out=100)
    for (i in seq(1, length(fitlist))) {
        if (fitlist[[i]]$errormodel == "yerrors") {
          plotwitherror(x = fitlist[[i]]$x, y = fitlist[[i]]$y,
                dy = fitlist[[i]]$dy, rep=TRUE, col=cols[i], ...)
        } else {
          plotwitherror(x = fitlist[[i]]$x, y = fitlist[[i]]$y,
                dy = fitlist[[i]]$dy, dx = fitlist[[i]]$dx, rep=TRUE, col=cols[i], ...)
        }
        
        lines(x=xseq, y=fitlist[[i]]$fn(par=fitlist[[i]]$t0, x=xseq, boot.R=0), col=cols[i])
        plotwitherror(x=0, y=fitlist[[i]]$t0[1], dy=fitlist[[i]]$se[1], col=cols[i], lwd=2, rep=T)
        if(drawlegend) do.call(what="legend", args=legendargs)
    }

}


minmax <- function(myvector) {
    lim <- c(min(myvector), max(myvector))
    return (lim*c(1-0.1*sign(lim[1]),1+0.1*sign(lim[2])))
}

errorpolygon <- function (X, fitresult, col.p, col.band = "gray",
                        polygon = TRUE, arlength = 0.1, pch = 1, ...) {
    # like plot of bootstrapfit, but with rep = TRUE
    prediction <- predict(fitresult, X)
    if (missing(pch)) {
        p.pch <- col.p
    } else {
        p.pch <- pch
    }
    if (fitresult$errormodel == "yerrors") {
      limits <- plotwitherror(x = fitresult$x, y = fitresult$y,
            dy = fitresult$dy, col = col.p, pch = p.pch, rep=TRUE, ...)
    } else {
      limits <- plotwitherror(x = fitresult$x, y = fitresult$y,
            dy = fitresult$dy, dx = fitresult$dx, col = col.p,
            pch = p.pch, rep=TRUE, ...)
    }
   ylim <- limits$ylim

    if (polygon) {
        polyval <- c(prediction$val + prediction$err,
                rev(prediction$val - prediction$err))
        if (any(polyval < ylim[1]) || any(polyval > ylim[2])) {
          polyval[polyval < ylim[1]] <- ylim[1]
          polyval[polyval > ylim[2]] <- ylim[2]
        }
        col.band <- "gray"
        pcol <- col2rgb(col.band, alpha = TRUE) / 255
        pcol[4] <- 0.65
        pcol <- rgb(red = pcol[1], green = pcol[2],
                    blue = pcol[3], alpha = pcol[4])
        polygon(x = c(X, rev(X)), y = polyval, col = pcol,
                lty = 0, lwd = 0.001, border = pcol)
    }
    lines(x = X, y = prediction$val, col = col.p, ...)
}
