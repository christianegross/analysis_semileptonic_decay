library("hadron")

cd <- read.table("/home/gross/Documents/heavymesons/scripts_christiane/extrapolate_sigma/complete_cB64_converted_new_sigmaepsilonto8amin1.csv")
cs <- read.table("/home/gross/Documents/heavymesons/scripts_christiane/extrapolate_sigma/cs_time_separationepsilonto8amin1e+0_converted.csv", header=TRUE)


cd <- cd[cd$th==6, ]
cs <- cs[cs$Nt==31, ]
#~ cd
#~ cs

#~ |Vcd|= 0.221±0.004 PDF chapter 12 Quark mixing matrix
#~ |Vcs|= 0.975±0.006

Vcd <- 0.221
Vcs <- 0.975


pdf("comparecscd.pdf", title="")


for (sys in c(T, F)) {
    title<- sprintf("B64, Nt31, th6, sys included %d", sys)
    maskcd <- cd$includesys==sys
    maskcs <- cs$includesys==sys
    plotwitherror(x=cd$iz[maskcd], y=cd$DGDq2[maskcd], dy=cd$dDGDq2[maskcd],
    col=1, pch=1, rep=FALSE,
    ylim = c(min(cd$DGDq2[maskcd] - cd$dDGDq2[maskcd], cs$DGDq2[maskcs] - cs$dDGDq2[maskcs]), max(cd$DGDq2[maskcd] + cd$dDGDq2[maskcd], cs$DGDq2[maskcs] + cs$dDGDq2[maskcs])),
    xlab="Z", ylab="24 pi^3 DGamma / Dq^2 [GeV^-3]", main=title)
    plotwitherror(x=cs$iz[maskcs], y=cs$DGDq2[maskcs], dy=cs$dDGDq2[maskcs],
    col=2, pch=2, rep=TRUE)
    legend(x="topleft", legend=c("cs", "cd"), pch=c(2, 1), col=c(2, 1))
}

for (sys in c(T, F)) {
    title<- sprintf("B64, Nt31, th6, sys included %d\nrescaled by |V|^2", sys)
    maskcd <- cd$includesys==sys
    maskcs <- cs$includesys==sys
    plotwitherror(x=cd$iz[maskcd], y=cd$DGDq2[maskcd]*Vcd^2, dy=cd$dDGDq2[maskcd]*Vcd^2,
    col=1, pch=1, rep=FALSE,
    ylim = c(min(cd$DGDq2[maskcd]*Vcd^2 - cd$dDGDq2[maskcd]*Vcd^2, cs$DGDq2[maskcs]*Vcs^2 - cs$dDGDq2[maskcs]*Vcs^2), 
            max(cd$DGDq2[maskcd]*Vcd^2 + cd$dDGDq2[maskcd]*Vcd^2, cs$DGDq2[maskcs]*Vcs^2 + cs$dDGDq2[maskcs]*Vcs^2)),
    xlab="Z", ylab="24 pi^3 DGamma / Dq^2 [GeV^-3]", main=title)
    plotwitherror(x=cs$iz[maskcs], y=cs$DGDq2[maskcs]*Vcs^2, dy=cs$dDGDq2[maskcs]*Vcs^2,
    col=2, pch=2, rep=TRUE)
    legend(x="topleft", legend=c("cs", "cd"), pch=c(2, 1), col=c(2, 1))
}


for (iz in c(3, 0, 1, 2)) {
    for (sys in c(T, F)) {
        title<- sprintf("B64, Nt31, th6, sys included %d, iZ %d", sys, iz)
        maskcd <- cd$iz==iz & cd$includesys==sys
        maskcs <- cs$iz==iz & cs$includesys==sys
        plotwitherror(x=cd$q[maskcd]^2, y=cd$DGDq2[maskcd], dy=cd$dDGDq2[maskcd],
        col=1, pch=1, rep=FALSE,
        ylim = c(min(cd$DGDq2[maskcd] - cd$dDGDq2[maskcd], cs$DGDq2[maskcs] - cs$dDGDq2[maskcs]), max(cd$DGDq2[maskcd] + cd$dDGDq2[maskcd], cs$DGDq2[maskcs] + cs$dDGDq2[maskcs])),
        xlab="q^2[GeV]", ylab="24 pi^3 DGamma / Dq^2 [GeV^-3]", main=title)
        plotwitherror(x=cs$q[maskcs]^2, y=cs$DGDq2[maskcs], dy=cs$dDGDq2[maskcs],
        col=2, pch=2, rep=TRUE)
        legend(x="topleft", legend=c("cs", "cd"), pch=c(2, 1), col=c(2, 1))
    }
}



for (iz in c(3, 0, 1, 2)) {
    for (sys in c(T, F)) {
        title<- sprintf("B64, Nt31, th6, sys included %d, iZ %d\nrescaled by |V|^2", sys, iz)
        maskcd <- cd$iz==iz & cd$includesys==sys
        maskcs <- cs$iz==iz & cs$includesys==sys
        plotwitherror(x=cd$q[maskcd]^2, y=cd$DGDq2[maskcd]*Vcd^2, dy=cd$dDGDq2[maskcd]*Vcd^2,
        col=1, pch=1, rep=FALSE,
        ylim = c(min(cd$DGDq2[maskcd]*Vcd^2 - cd$dDGDq2[maskcd]*Vcd^2, cs$DGDq2[maskcs]*Vcs^2 - cs$dDGDq2[maskcs]*Vcs^2), 
                max(cd$DGDq2[maskcd]*Vcd^2 + cd$dDGDq2[maskcd]*Vcd^2, cs$DGDq2[maskcs]*Vcs^2 + cs$dDGDq2[maskcs]*Vcs^2)),
        xlab="q^2[GeV]", ylab="24 pi^3 DGamma / Dq^2 [GeV^-3]", main=title)
        plotwitherror(x=cs$q[maskcs]^2, y=cs$DGDq2[maskcs]*Vcs^2, dy=cs$dDGDq2[maskcs]*Vcs^2,
        col=2, pch=2, rep=TRUE)
        legend(x="topleft", legend=c("cs", "cd"), pch=c(2, 1), col=c(2, 1))
    }
}

