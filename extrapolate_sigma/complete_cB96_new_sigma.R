library("hadron")
source("/home/gross/Documents/heavymesons/scripts_christiane/extrapolate_sigma_templates/calc_DGDq2.R")

#~ data <- read_in_DGDq2(filename="/home/gross/Documents/heavymesons/data/cB211ab.07.64/t48_64/outputDGammaDq2/DGammaDq2.bin", write=FALSE)
#~ data <- read_in_DGDq2(filename="/home/gross/Documents/heavymesons/data/cB211ab.07.64/th6_t48/outputDGammaDq2_10sigma/DGammaDq2.bin", write=FALSE)

## for different possible combinations of covariance diagonal or full and error statistical or total, 
## we read in the results of the analysis for different parameters (time, set (Z0, Z1, Z2, full) and comb (full, vpar, apar, vperp, aperp).
## for each parameter combination, we extract the epsilons used, and the results for DG from the data, and fit it with a parabolic function to get the result for epsilon=0.
## We can choose to do the fit using the cavariance matrix estimated from the bootstrap samples or only from the errors.


comment <- "epsilonto8amin1"

if(TRUE) {
fnqar <- function(par, x, boot.R, ...) par[1] + par[2] * x**2 + par[3] * x**4
fnparnolin <- function(par, x, boot.R, ...) par[1] + par[2] * x**2
fnpar <- function(par, x, boot.R, ...) par[1] + par[2] * x + par[3] * x**2
fnlin <- function(par, x, boot.R, ...) par[1] + par[2] * x
fncon <- function(par, x, boot.R, ...) par[1]
fneight <- function(par, x, boot.R, ...) par[1] + par[2] * x**2 + par[3] * x**4 + par[4] * x**6 + par[5] * x**8

fnlist <- list(fncon, fnlin, fnpar, fnparnolin, fnqar)



pdf(sprintf("complete_cB96_new_sigma%s.pdf", comment), title="")

dosys <- T
uplim <- 18
lowlim <- 3
degree <- 4
## largest permissible A/A0_ref
maxAmin <- 1


files <- c("th1/new_sigma/outputDGammaDq2/DGammaDq2.bin",
           "th2/new_sigma/outputDGammaDq2/DGammaDq2.bin",
           "th3/new_sigma/outputDGammaDq2/DGammaDq2.bin",
           "th4/new_sigma/outputDGammaDq2/DGammaDq2.bin",
           "th5/new_sigma/outputDGammaDq2/DGammaDq2.bin",
           "th6/new_sigma/outputDGammaDq2/DGammaDq2.bin",
           "th7/new_sigma/outputDGammaDq2/DGammaDq2.bin",
           "th8/new_sigma/outputDGammaDq2/DGammaDq2.bin",
           "th9/new_sigma/outputDGammaDq2/DGammaDq2.bin",
           "th9.5/new_sigma/outputDGammaDq2/DGammaDq2.bin")

tsnk <- c(56, 56, 56, 56, 56, 56, 56, 56, 56, 56)
Nt <- c(31, 31, 31, 31, 31, 31, 31, 31, 31, 31)
th <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 9.5)
nerr <- c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2)

savename <- sprintf("complete_cB96_new_sigma%s", comment)
par.guess <- rep(1, 5)

determineDGDq2(resultpath = "/home/gross/Documents/heavymesons/data/cB211.07.96/", filenames = files, 
        tsnk = tsnk, Nt = Nt, th = th, nerr = nerr, amin = maxAmin, savename = savename,
        fitfn = fneight, par.guess = par.guess)

}

if(TRUE) {
a <- 0.07957 / 0.1973269804 # fm / hbarc = GeV^-1
am_H <- 7.996044e-01 # dimensionless
am_H <- 7.996939e-01 # dimensionless

conversionfactor <- (am_H/a)^3/2/pi

print(conversionfactor)

savename <- sprintf("complete_cB96_new_sigma%s.csv", comment)
result <- read.table(savename, header=TRUE)
pdf(sprintf("complete_cB96_DqDq2_new_sigma%s.pdf", comment), title="")

result[, c("DGDq2", "dDGDq2")] <- result[, c("DGDq2", "dDGDq2")] * conversionfactor
result$q <- result$w * am_H / a
result$agev <- rep(a, length(result$q))

mask <- result$iz==3 & result$icomb==0 & result$includesys==TRUE
print(result[mask, ])

plotwitherror(x=result$q[mask & result$nerr==2]^2, y=result$DGDq2[mask & result$nerr==2], dy=result$dDGDq2[mask & result$nerr==2],
main="differential decay rate for tsink=56, tins=44, including systematical error", xlab="q^2 [GeV^2]", ylab="24 pi^3 DGamma / Dq^2 [GeV^-3]",
ylim=c(0, max(result$DGDq2[mask] + result$dDGDq2[mask])))
legend(x="topleft", legend=c("nerr 2"), col=c(1), pch=c(1))

plotwitherror(x=result$q[mask & result$nerr==2]^2, y=result$dDGDq2[mask & result$nerr==2],
main="errors differential decay rate for tsink=56, tins=44, including systematical error", xlab="q^2 [GeV^2]", ylab="24 pi^3 d(DGamma / Dq^2) [GeV^-3]")
plotwitherror(x=result$q[mask & result$nerr==1]^2, y=result$dDGDq2[mask & result$nerr==1],
rep=TRUE, col=2, pch=2)
legend(x="topleft", legend=c("nerr 2", "nerr 1"), col=c(1, 2), pch=c(1, 2))

plotwitherror(x=result$q[mask & result$nerr==2 & result$dDGDq2 < 1]^2, y=result$dDGDq2[mask & result$nerr==2 & result$dDGDq2 < 1],
main="errors differential decay rate for tsink=56, tins=44, including systematical error", xlab="q^2 [GeV^2]", ylab="24 pi^3 d(DGamma / Dq^2) [GeV^-3]")
legend(x="topleft", legend=c("nerr 2"), col=c(1), pch=c(1))

plotwitherror(x=result$q[mask & result$nerr==2]^2, y=result$dDGDq2[mask & result$nerr==2] / abs(result$DGDq2[mask & result$nerr==2]),
main="rel. errors differential decay rate for tsink=56, tins=44, including systematical error", xlab="q^2 [GeV^2]", ylab="24 pi^3 d(DGamma / Dq^2) [GeV^-3]")
legend(x="topleft", legend=c("nerr 2"), col=c(1), pch=c(1))

mask <- result$iz==3 & result$icomb==0 & result$includesys==FALSE
print(result[mask, ])

print(max(result$DGDq2[mask] + result$dDGDq2[mask]))

plotwitherror(x=result$q[mask & result$nerr==2]^2, y=result$DGDq2[mask & result$nerr==2], dy=result$dDGDq2[mask & result$nerr==2],
main="differential decay rate for tsink=56, tins=44, excluding systematical error", xlab="q^2 [GeV^2]", ylab="24 pi^3 DGamma / Dq^2 [GeV^-3]",
ylim=c(0, max(result$DGDq2[mask] + result$dDGDq2[mask])))
legend(x="topleft", legend=c("nerr 2"), col=c(1), pch=c(1))

plotwitherror(x=result$q[mask & result$nerr==2 & result$dDGDq2 < 1]^2, y=result$dDGDq2[mask & result$nerr==2 & result$dDGDq2 < 1],
main="errors differential decay rate for tsink=56, tins=44, excluding systematical error", xlab="q^2 [GeV^2]", ylab="24 pi^3 d(DGamma / Dq^2) [GeV^-3]")
legend(x="topleft", legend=c("nerr 2"), col=c(1), pch=c(1))

plotwitherror(x=result$q[mask & result$nerr==2]^2, y=result$dDGDq2[mask & result$nerr==2],
main="errors differential decay rate for tsink=56, tins=44, excluding systematical error", xlab="q^2 [GeV^2]", ylab="24 pi^3 d(DGamma / Dq^2) [GeV^-3]")
legend(x="topleft", legend=c("nerr 2"), col=c(1), pch=c(1))

plotwitherror(x=result$q[mask & result$nerr==2]^2, y=result$dDGDq2[mask & result$nerr==2] / abs(result$DGDq2[mask & result$nerr==2]),
main="rel. errors differential decay rate for tsink=56, tins=44, excluding systematical error", xlab="q^2 [GeV^2]", ylab="24 pi^3 d(DGamma / Dq^2) [GeV^-3]")
legend(x="topleft", legend=c("nerr 2"), col=c(1), pch=c(1))


## relative errors
## autocorrelation times

write.table(result, sprintf("complete_cB96_converted_new_sigma%s.csv", comment))

}

## Extract minimum A/A0 for each Z, check if it is permissible. Save this in array of dimensions #files*#epsilon, by using AND
## Only then do loop over Z3
