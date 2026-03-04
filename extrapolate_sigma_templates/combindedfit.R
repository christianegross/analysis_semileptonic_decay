library("hadron")

predict.bootstrapfit.comb <- function (object, x, error = object$error.function, maskfn, ...) {
  ## to include additional parameter to x$fn originally given as ... to
  ## bootstrap.nlsfit requires some pull-ups
  npar <- length(object$par.guess)
  val <- do.call(object$fn, c(list(par = object$t0[1:npar], x = x, boot.r = 0, maskfn=maskfn)))
  
  prediction <- list(x = x, val = val)
  
  if(!is.null(object$t)) {
    ## error band
    ## define a dummy function to be used in apply
    prediction_boot_fn <- function (boot.r) {
      par <- object$t[boot.r, 1:npar, drop = FALSE]
      do.call(object$fn, c(list(par = par, x = x, boot.r = boot.r, maskfn=maskfn)))
    }
    prediction_boot <- do.call(rbind, lapply(1:nrow(object$t), prediction_boot_fn))
    prediction$boot <- prediction_boot
    
    err <- apply(prediction_boot, 2, error, na.rm = TRUE)
    stopifnot(length(err) == length(x))
    prediction$err <- err
  }
  
  return (invisible(prediction))
}
# pdf("trycombined.pdf", title="")

plotcombined <- function(fitresult, cols=c("black", "red"), pchs=c(1, 2), xlim, ylim, polygon=TRUE, ...) {
if(missing(xlim)) xlim <- range(fitresult$x)
if(missing(ylim)) ylim <- range(na.omit(fitresult$y+fitresult$dy, fitresult$y-fitresult$dy))
xseq <- seq(from=xlim[1]*(1-2*sign(xlim[1])), to=xlim[2]*(1+0.5*sign(xlim[2])), length.out=300)
y1 <- predict.bootstrapfit.comb(fitresult, x=xseq, maskfn=T)
y2 <- predict.bootstrapfit.comb(fitresult, x=xseq, maskfn=F)

plotwitherror(x=fitresult$x, y=fitresult$y, dy=fitresult$dy, col=ifelse(fitresult$tofn$maskfn, cols[1], cols[2]), pch=ifelse(fitresult$tofn$maskfn, pchs[1], pchs[2]), xlim=xlim, ylim=ylim, ...)
points(x=xseq, y=y1$val, type = "l", col=cols[1])
points(x=xseq, y=y2$val, type = "l", col=cols[2])

if(polygon) {
polyval1 <- c(y1$val + y1$err, rev(y1$val - y1$err))
polyval2 <- c(y2$val + y2$err, rev(y2$val - y2$err))

pcol1 <- col2rgb(cols[1], alpha=TRUE)/255 
pcol1[4] <- 0.2
pcol1 <- rgb(red=pcol1[1],green=pcol1[2],blue=pcol1[3],alpha=pcol1[4])
pcol2 <- col2rgb(cols[2], alpha=TRUE)/255 
pcol2[4] <- 0.2
pcol2 <- rgb(red=pcol2[1],green=pcol2[2],blue=pcol2[3],alpha=pcol2[4])
polygon(x=c(xseq, rev(xseq)), y=polyval1, col=pcol1, lty=0, lwd=0.001, border=pcol1)
polygon(x=c(xseq, rev(xseq)), y=polyval2, col=pcol2, lty=0, lwd=0.001, border=pcol2)
}
plotwitherror(x=fitresult$x, y=fitresult$y, dy=fitresult$dy, col=ifelse(fitresult$tofn$maskfn, cols[1], cols[2]), pch=ifelse(fitresult$tofn$maskfn, pchs[1], pchs[2]), xlim=xlim, ylim=ylim, rep=T, ...)

}


#~ x <- c(1:10, 1:10)
#~ y1 <- 1:10*2+rnorm(10, mean=0, sd=0.5)
#~ y2 <- 1:10*3+rnorm(10, mean=0, sd=0.1)
#~ y <- c(y1, y2)
#~ bsamples <- parametric.bootstrap(100, x=c(y1, y2), dx=rep(0.5, 20), seed=1234)
#~ fn <- function(par, x, boot.r, maskfn, ...) {
#~   par[1] + par[2]*x*maskfn + par[3]*x*(!maskfn)
#~ }
#~ fit.result <- bootstrap.nlsfit(fn=fn, par.guess=c(1, 1, 1), y=y, x=x, bsamples = bsamples, maskfn=c(rep(T, 10), rep(F, 10)))
#~ summary(fit.result)
#~ plotcombined(fit.result, cols=c("black", "orange"), pch=c(21, 21), xlim=c(0, 10), ylim=c(0, 30))
#~ plotwitherror(fit.result$x, fit.result$y, fit.result$dy, rep=T, col=ifelse(fit.result$tofn$maskfn, "black", "orange"))
