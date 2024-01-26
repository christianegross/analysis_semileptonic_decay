library("hadron")

B64 <- readRDS("complete_cB64_new_sigmaepsilonto8amin1.RDS")
B96 <- readRDS("complete_cB96_new_sigmaepsilonto8amin1.RDS")
C80 <- readRDS("complete_cC80_new_sigmaepsilonto8amin1.RDS")
D96 <- readRDS("complete_cD96_new_sigmaepsilonto8amin1.RDS")

B64table <- read.table("complete_cB64_converted_new_sigmaepsilonto8amin1.csv")
B96table <- read.table("complete_cB96_converted_new_sigmaepsilonto8amin1.csv")
C80table <- read.table("complete_cC80_converted_new_sigmaepsilonto8amin1.csv")
D96table <- read.table("complete_cD96_converted_new_sigmaepsilonto8amin1.csv")

B64$L <- rep(64, length(B64$Nt))
B96$L <- rep(96, length(B96$Nt))
C80$L <- rep(80, length(C80$Nt))
D96$L <- rep(96, length(D96$Nt))

bsamples <- array(NA, dim=c(1000, 4))
y <- c()
dy <- c()
x <- c()

enslist <- list(B64, B96, C80, D96)
tablelist <- list(B64table, B96table, C80table, D96table)


fnlin <- function(par, x, boot.R, ...) par[1] + par[2] * x
fncon <- function(par, x, boot.R, ...) par[1]

fitfn <- fnlin
par.guess <- c(1, 1)
comment <- "linear"

res <- data.frame(th=NA, sys=NA, iz=NA, lim=NA, dlim=NA, q=NA)


pcol <- col2rgb("gray", alpha=TRUE)/255 
pcol[4] <- 0.65
pcol <- rgb(red=pcol[1],green=pcol[2],blue=pcol[3],alpha=pcol[4])

pdf(sprintf("contlimit%s.pdf", comment))

for (th in c(1, 2, 3, 4, 5, 6, 7, 8, 9, 9.5)){
    print(sprintf("%f %f %f %f", B64table$q[abs(B64table$th - th) < 1e-2][1], B96table$q[abs(B96table$th - th) < 1e-2][1], C80table$q[abs(C80table$th - th) < 1e-2][1], D96table$q[abs(D96table$th - th) < 1e-2][1]))
    print(sprintf("%e %e %e", B64table$q[abs(B64table$th - th) < 1e-2][1] - B96table$q[abs(B96table$th - th) < 1e-2][1], 
                              B64table$q[abs(B64table$th - th) < 1e-2][1] - C80table$q[abs(C80table$th - th) < 1e-2][1], 
                              B64table$q[abs(B64table$th - th) < 1e-2][1] - D96table$q[abs(D96table$th - th) < 1e-2][1]))
    for (sys in c(T, F)) {
        for (iz in c(3, 0, 1, 2)) {
            title <- sprintf("th%s q^2 %.3f GeV^2 iz %d include sys %d", th, B64table$q[abs(B64table$th - th) < 1e-2][1]^2, iz, sys)
#~             print(title)
            for (index in seq(1, 4)) {
                mask <- abs(enslist[[index]]$th - th) < 1e-2 & enslist[[index]]$iz == iz & enslist[[index]]$includesys==sys
                masktable <- abs(tablelist[[index]]$th - th) < 1e-2 & tablelist[[index]]$iz == iz & tablelist[[index]]$includesys==sys
                y[index] <- enslist[[index]]$DGDq2[mask]
                dy[index] <- enslist[[index]]$dDGDq2[mask]
                bsamples[, index] <- enslist[[index]]$dat[, mask]
                x[index] <- tablelist[[index]]$agev[masktable]
            }
            
            # comvert x from gev to gm 
            x <- x * 0.1973269804
            
            fitresult <- bootstrap.nlsfit(x=x^2, y=y, dy=dy, bsamples=bsamples, fn=fitfn, par.guess=par.guess)
            
            if(!inherits(fitresult, "try-error")) {
                if(length(par.guess)==1) {
                    plotwitherror(x=x^2, y=y, dy=dy, xlab="a[fm]^2", ylab="24 pi^3 DGamma / Dq^2 [GeV^-3]", main=title, xlim=c(-0.001, 0.007))
                    polygon(x=c(-0.002, 0.01, 0.01, -0.002), y=c(rep(fitresult$t0[1] - fitresult$se[1], 2), rep(fitresult$t0[1] + fitresult$se[1], 2)), 
                    col=pcol)
                    lines(x=c(-0.002, 0.01), y=rep(fitresult$t0[1], 2))
                    plotwitherror(x=x^2, y=y, dy=dy, rep=TRUE)
                } else {
                    plot(fitresult, xlab="a[fm]^2", ylab="24 pi^3 DGamma / Dq^2 [GeV^-3]", main=title, xlim=c(-0.001, 0.007), plot.range=2*c(-0.001, 0.007))
                }
                plotwitherror(x=0, y=fitresult$t0[1], dy=fitresult$se[1], rep=TRUE, col="red")
                res <- rbind(res, data.frame(th=th, iz=iz, sys=sys, lim=fitresult$t0[1], dlim=fitresult$se[1], q=B64table$q[abs(B64table$th - th) < 1e-2][1]))
            } else {
                plotwitherror(x=x, y=y, dy=dy, xlab="a[fm]^2", ylab="24 pi^3 DGamma / Dq^2 [GeV^-3]", main=title)
            }
        }
    }
}

write.table(res, sprintf("contlimit%s.csv", comment))

for (sys in c(T, F)) {
    for (iz in c(3, 0, 1, 2)) {
        mask <- res$iz==iz & res$sys==sys
        plotwitherror(x=res$q[mask]^2, y=res$lim[mask], dy=res$dlim[mask], 
        xlab="q^2 [GeV^2]", ylab="24 pi^3 DGamma / Dq^2 [GeV^-3]", 
        main=paste("continuum limit sys included ", sys, "iZ", iz))
    }
}
