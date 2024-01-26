library("hadron")
pdf("comparevolume_new_sigma_epsilonto8_amin1.pdf", title="")


cB64 <- read.table("complete_cB64_converted_new_sigmaepsilonto8amin1.csv")
cB96 <- read.table("complete_cB96_converted_new_sigmaepsilonto8amin1.csv")
cC80 <- read.table("complete_cC80_converted_new_sigmaepsilonto8amin1.csv")
cD96 <- read.table("complete_cD96_converted_new_sigmaepsilonto8amin1.csv")

cB64 <- na.omit(cB64)
cB96 <- na.omit(cB96)
cC80 <- na.omit(cC80)
cD96 <- na.omit(cD96)

cB64 <- cB64[order(cB64$iz, cB64$includesys, cB64$icomb, cB64$q), ]
cB96 <- cB96[order(cB96$iz, cB96$includesys, cB96$icomb, cB96$q), ]
cC80 <- cC80[order(cC80$iz, cC80$includesys, cC80$icomb, cC80$q), ]
cD96 <- cD96[order(cD96$iz, cD96$includesys, cD96$icomb, cD96$q), ]

sumup <- TRUE

if(sumup) {
sumres <- data.frame(q=rep(0, 20), sys=c(rep(TRUE, 10), rep(FALSE, 10)), B64DGDq2=rep(0, 20), B64dDGDq2=rep(0, 20), 
            B96DGDq2=rep(0, 20), B96dDGDq2=rep(0, 20), C80DGDq2=rep(0, 20), 
            C80dDGDq2=rep(0, 20), D96DGDq2=rep(0, 20), D96dDGDq2=rep(0, 20))
}

for(iz in c(3, 0, 1, 2)) {

maskB64 <- cB64$iz==iz & cB64$icomb==0 & cB64$includesys==TRUE
maskB96 <- cB96$iz==iz & cB96$icomb==0 & cB96$includesys==TRUE
maskC80 <- cC80$iz==iz & cC80$icomb==0 & cC80$includesys==TRUE
maskD96 <- cD96$iz==iz & cD96$icomb==0 & cD96$includesys==TRUE

print(cB64$q[maskB64])

plotwitherror(x=cB64$q[maskB64]^2, y=cB64$DGDq2[maskB64], dy=cB64$dDGDq2[maskB64],
main=paste("diff. decay rate for tsink=56, tins=44, sys+stat, z =", iz), xlab="q^2 [GeV^2]", ylab="24 pi^3 DGamma / Dq^2 [GeV^-3]",
ylim=c(0, max(cB64$DGDq2[maskB64] + cB64$dDGDq2[maskB64], cB96$DGDq2[maskB96] + cB96$dDGDq2[maskB96])))
plotwitherror(x=cB96$q[maskB96]^2, y=cB96$DGDq2[maskB96], dy=cB96$dDGDq2[maskB96],
rep=TRUE, col=2, pch=2)
plotwitherror(x=cC80$q[maskC80]^2, y=cC80$DGDq2[maskC80], dy=cC80$dDGDq2[maskC80],
rep=TRUE, col=5, pch=5)
plotwitherror(x=cD96$q[maskD96]^2, y=cD96$DGDq2[maskD96], dy=cD96$dDGDq2[maskD96],
rep=TRUE, col=4, pch=4)
legend(x="topleft", legend=c("cB64", "cB96", "cC80", "cD96"), col=c(1, 2, 5, 4), pch=c(1, 2, 5, 4))


if(sumup) {
if(iz==0) { sumres$q[sumres$sys==TRUE] <- cB64$q[maskB64] }
if(iz != 3) {
    sumres$B64DGDq2[sumres$sys==TRUE]  <- sumres$B64DGDq2[sumres$sys==TRUE]  + cB64$DGDq2[maskB64]
    sumres$B64dDGDq2[sumres$sys==TRUE] <- sumres$B64dDGDq2[sumres$sys==TRUE] + cB64$dDGDq2[maskB64]^2
    sumres$B96DGDq2[sumres$sys==TRUE]  <- sumres$B96DGDq2[sumres$sys==TRUE]  + cB96$DGDq2[maskB96]
    sumres$B96dDGDq2[sumres$sys==TRUE] <- sumres$B96dDGDq2[sumres$sys==TRUE] + cB96$dDGDq2[maskB96]^2
    sumres$C80DGDq2[sumres$sys==TRUE]  <- sumres$C80DGDq2[sumres$sys==TRUE]  + cC80$DGDq2[maskC80]
    sumres$C80dDGDq2[sumres$sys==TRUE] <- sumres$C80dDGDq2[sumres$sys==TRUE] + cC80$dDGDq2[maskC80]^2
    sumres$D96DGDq2[sumres$sys==TRUE]  <- sumres$D96DGDq2[sumres$sys==TRUE]  + cD96$DGDq2[maskD96]
    sumres$D96dDGDq2[sumres$sys==TRUE] <- sumres$D96dDGDq2[sumres$sys==TRUE] + cD96$dDGDq2[maskD96]^2
}}

maskB64 <- cB64$iz==iz & cB64$icomb==0 & cB64$includesys==FALSE
maskB96 <- cB96$iz==iz & cB96$icomb==0 & cB96$includesys==FALSE
maskC80 <- cC80$iz==iz & cC80$icomb==0 & cC80$includesys==FALSE
maskD96 <- cD96$iz==iz & cD96$icomb==0 & cD96$includesys==FALSE

plotwitherror(x=cB64$q[maskB64]^2, y=cB64$DGDq2[maskB64], dy=cB64$dDGDq2[maskB64],
main=paste("diff. decay rate for tsink=56, tins=44, stat, z =", iz), xlab="q^2 [GeV^2]", ylab="24 pi^3 DGamma / Dq^2 [GeV^-3]",
ylim=c(0, max(cB64$DGDq2[maskB64] + cB64$dDGDq2[maskB64], cB96$DGDq2[maskB96] + cB96$dDGDq2[maskB96])))
plotwitherror(x=cB96$q[maskB96]^2, y=cB96$DGDq2[maskB96], dy=cB96$dDGDq2[maskB96],
rep=TRUE, col=2, pch=2)
plotwitherror(x=cC80$q[maskC80]^2, y=cC80$DGDq2[maskC80], dy=cC80$dDGDq2[maskC80],
rep=TRUE, col=5, pch=5)
plotwitherror(x=cD96$q[maskD96]^2, y=cD96$DGDq2[maskD96], dy=cD96$dDGDq2[maskD96],
rep=TRUE, col=4, pch=4)
legend(x="topleft", legend=c("cB64", "cB96", "cC80", "cD96"), col=c(1, 2, 5, 4), pch=c(1, 2, 5, 4))



if(sumup) {
if(iz==0) { sumres$q[sumres$sys==FALSE] <- cB64$q[maskB64] }
if(iz != 3) {
    sumres$B64DGDq2[sumres$sys==FALSE]  <- sumres$B64DGDq2[sumres$sys==FALSE]  + cB64$DGDq2[maskB64]
    sumres$B64dDGDq2[sumres$sys==FALSE] <- sumres$B64dDGDq2[sumres$sys==FALSE] + cB64$dDGDq2[maskB64]^2
    sumres$B96DGDq2[sumres$sys==FALSE]  <- sumres$B96DGDq2[sumres$sys==FALSE]  + cB96$DGDq2[maskB96]
    sumres$B96dDGDq2[sumres$sys==FALSE] <- sumres$B96dDGDq2[sumres$sys==FALSE] + cB96$dDGDq2[maskB96]^2
    sumres$C80DGDq2[sumres$sys==FALSE]  <- sumres$C80DGDq2[sumres$sys==FALSE]  + cC80$DGDq2[maskC80]
    sumres$C80dDGDq2[sumres$sys==FALSE] <- sumres$C80dDGDq2[sumres$sys==FALSE] + cC80$dDGDq2[maskC80]^2
    sumres$D96DGDq2[sumres$sys==FALSE]  <- sumres$D96DGDq2[sumres$sys==FALSE]  + cD96$DGDq2[maskD96]
    sumres$D96dDGDq2[sumres$sys==FALSE] <- sumres$D96dDGDq2[sumres$sys==FALSE] + cD96$dDGDq2[maskD96]^2
}}


}

## sum of Z0, Z1, Z2


if(sumup) {
print(sumres)

print(sumres$q[sumres$sys==TRUE]^2)
print(sumres$B64DGDq2[sumres$sys==TRUE])
print(sqrt(sumres$B64dDGDq2[sumres$sys==TRUE]))

plotwitherror(x=sumres$q[sumres$sys==TRUE]^2, y=sumres$B64DGDq2[sumres$sys==TRUE], dy=sqrt(sumres$B64dDGDq2[sumres$sys==TRUE]),
main=paste("diff. decay rate sys+stat, Z0+Z1+Z2"), xlab="q^2 [GeV^2]", ylab="24 pi^3 DGamma / Dq^2 [GeV^-3]")
plotwitherror(x=sumres$q[sumres$sys==TRUE]^2, y=sumres$B96DGDq2[sumres$sys==TRUE], dy=sqrt(sumres$B96dDGDq2[sumres$sys==TRUE]),
rep=TRUE, col=2, pch=2)
plotwitherror(x=sumres$q[sumres$sys==TRUE]^2, y=sumres$C80DGDq2[sumres$sys==TRUE], dy=sqrt(sumres$C80dDGDq2[sumres$sys==TRUE]),
rep=TRUE, col=5, pch=5)
plotwitherror(x=sumres$q[sumres$sys==TRUE]^2, y=sumres$D96DGDq2[sumres$sys==TRUE], dy=sqrt(sumres$D96dDGDq2[sumres$sys==TRUE]),
rep=TRUE, col=4, pch=4)
legend(x="topleft", legend=c("cB64", "cB96", "cC80", "cD96"), col=c(1, 2, 5, 4), pch=c(1, 2, 5, 4))

plotwitherror(x=sumres$q[sumres$sys==FALSE]^2, y=sumres$B64DGDq2[sumres$sys==FALSE], dy=sqrt(sumres$B64dDGDq2[sumres$sys==FALSE]),
main=paste("diff. decay rate stat, Z0+Z1+Z2"), xlab="q^2 [GeV^2]", ylab="24 pi^3 DGamma / Dq^2 [GeV^-3]")
plotwitherror(x=sumres$q[sumres$sys==FALSE]^2, y=sumres$B96DGDq2[sumres$sys==FALSE], dy=sqrt(sumres$B96dDGDq2[sumres$sys==FALSE]),
rep=TRUE, col=2, pch=2)
plotwitherror(x=sumres$q[sumres$sys==FALSE]^2, y=sumres$C80DGDq2[sumres$sys==FALSE], dy=sqrt(sumres$C80dDGDq2[sumres$sys==FALSE]),
rep=TRUE, col=5, pch=5)
plotwitherror(x=sumres$q[sumres$sys==FALSE]^2, y=sumres$D96DGDq2[sumres$sys==FALSE], dy=sqrt(sumres$D96dDGDq2[sumres$sys==FALSE]),
rep=TRUE, col=4, pch=4)
legend(x="topleft", legend=c("cB64", "cB96", "cC80", "cD96"), col=c(1, 2, 5, 4), pch=c(1, 2, 5, 4))

}
