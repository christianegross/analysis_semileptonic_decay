library("hadron")
source("/home/gross/Documents/heavymesons/scripts_christiane/extrapolate_sigma_templates/calc_DGDq2.R")

respath <- "/home/gross/Documents/heavymesons/data/cs/cB211.07.64/test_time_separation"

files <- c("th6_t56_44/Nt31/outputDGammaDq2/DGammaDq2.bin",
           "th6_t56_44/Nt30/outputDGammaDq2/DGammaDq2.bin",
           "th6_t56_44/Nt29/outputDGammaDq2/DGammaDq2.bin",
           "th6_t56_44/Nt28/outputDGammaDq2/DGammaDq2.bin",
           "th6_t56_44/Nt27/outputDGammaDq2/DGammaDq2.bin",
           "th6_t56_44/Nt26/outputDGammaDq2/DGammaDq2.bin",
           "th6_t56_44/Nt25/outputDGammaDq2/DGammaDq2.bin",
           "th6_t52_40/Nt27/outputDGammaDq2/DGammaDq2.bin",
           "th6_t52_40/Nt26/outputDGammaDq2/DGammaDq2.bin",
           "th6_t52_40/Nt25/outputDGammaDq2/DGammaDq2.bin",
           "th6_t52_40/Nt24/outputDGammaDq2/DGammaDq2.bin",
           "th6_t52_40/Nt23/outputDGammaDq2/DGammaDq2.bin",
           "th6_t52_40/Nt22/outputDGammaDq2/DGammaDq2.bin"
           )
           
Nt <- c(seq(31, 25), seq(27, 22))
tsnk <- c(rep(56, 7), rep(52, 6))
th <- rep(6, 13)
nerr <- rep(2, 13)

amin <- 1e+0
comment <- "epsilonto8amin1e+0"
savename <- sprintf("cs_time_separation%s", comment)

if(TRUE) {
pdf(sprintf("cs_time_separation%s.pdf", comment), title="")

fnqar <- function(par, x, boot.R, ...) par[1] + par[2] * x**2 + par[3] * x**4
par.guess <- rep(1, 3)
fnsix <- function(par, x, boot.R, ...) par[1] + par[2] * x**2 + par[3] * x**4 + par[4] * x**6
par.guess <- rep(1, 4)

fneight <- function(par, x, boot.R, ...) par[1] + par[2] * x**2 + par[3] * x**4 + par[4] * x**6 + par[5] * x**8
par.guess <- rep(1, 5)

determineDGDq2(resultpath=respath, filenames=files, tsnk=tsnk, Nt=Nt, 
        th=th, nerr=nerr, amin=amin, fitfn=fneight, par.guess=par.guess,
        savename=savename, epsuplim=0)

#~ ## 
#~ determineDGDq2 <- function(resultpath, filenames, tsnk, Nt, th, 
#~         nerr, amin, savename, fitfn, par.guess, isets=c(0), icomb=c(0)) 
}


pdf(sprintf("cs_time_separation_DGDq2%s.pdf", comment), title="")


result <- read.table(paste0(savename, ".csv"), header=TRUE)
names(result)


a <- 0.07957 / 0.1973269804 # fm / hbarc = GeV^-1
am_H <- 7.996044e-01 # dimensionless
conversionfactor <- (am_H/a)^3/2/pi
result[, c("DGDq2", "dDGDq2")] <- result[, c("DGDq2", "dDGDq2")] * conversionfactor
result$q <- result$w * am_H / a
result$agev <- rep(a, length(result$q))


for(iz in c(3, 0, 1, 2)) {
    mask <- result$iz==3 & result$icomb==0 & result$includesys==TRUE
    
    
    plotwitherror(x=result$Nt[result$tsnk==56 & mask], y=result$DGDq2[result$tsnk==56 & mask], dy=result$dDGDq2[result$tsnk==56 & mask],
    xlab="Nt", ylab="DGDq2", main=paste("B64 th6, sys+stat, z =", iz), xlim=c(min(result$Nt), max(result$Nt)))
    
    plotwitherror(x=result$Nt[result$tsnk==52 & mask]+0.1, y=result$DGDq2[result$tsnk==52 & mask], dy=result$dDGDq2[result$tsnk==52 & mask],
    rep=TRUE, col=2, pch=2)
    legend(x="top", legend=c("tsink=56", "tsink=52"), pch=c(1, 2), col=c(1, 2))
    
    mask <- result$iz==3 & result$icomb==0 & result$includesys==FALSE
    
    plotwitherror(x=result$Nt[result$tsnk==56 & mask], y=result$DGDq2[result$tsnk==56 & mask], dy=result$dDGDq2[result$tsnk==56 & mask],
    xlab="Nt", ylab="DGDq2", main=paste("B64 th6, stat, z =", iz), xlim=c(min(result$Nt), max(result$Nt)))
    
    plotwitherror(x=result$Nt[result$tsnk==52 & mask]+0.1, y=result$DGDq2[result$tsnk==52 & mask], dy=result$dDGDq2[result$tsnk==52 & mask],
    rep=TRUE, col=2, pch=2)
    legend(x="top", legend=c("tsink=56", "tsink=52"), pch=c(1, 2), col=c(1, 2))
}


write.table(result, sprintf("%s_converted.csv", savename))

