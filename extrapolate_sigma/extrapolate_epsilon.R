library("hadron")
source("read_in_binary.R")

#~ data <- read_in_DGDq2(filename="/home/gross/Documents/heavymesons/data/cB211ab.07.64/t48_64/outputDGammaDq2/DGammaDq2.bin", write=FALSE)
#~ data <- read_in_DGDq2(filename="/home/gross/Documents/heavymesons/data/cB211ab.07.64/th6_t48/outputDGammaDq2_10sigma/DGammaDq2.bin", write=FALSE)

## for different possible combinations of covariance diagonal or full and error statistical or total, 
## we read in the results of the analysis for different parameters (time, set (Z0, Z1, Z2, full) and comb (full, vpar, apar, vperp, aperp).
## for each parameter combination, we extract the epsilons used, and the results for DG from the data, and fit it with a parabolic function to get the result for epsilon=0.
## We can choose to do the fit using the cavariance matrix estimated from the bootstrap samples or only from the errors.

if(TRUE) {
fnpar <- function(par, x, boot.R, ...) par[1] + par[2] * x + par[3] * x**2
fnlin <- function(par, x, boot.R, ...) par[1] + par[2] * x
fncon <- function(par, x, boot.R, ...) par[1]


covfull <- T
fitcorr <- F

if (!covfull) ptitle <- "extrapolate_epsilon_covdiag"
if (covfull)  ptitle <- "extrapolate_epsilon_covfull"
if (fitcorr) ptitle <- paste0(ptitle, "fitcorr.pdf")
if (!fitcorr) ptitle <- paste0(ptitle, "fituncr.pdf")

pdf(ptitle, title="")

dosys <- T
uplim <- 10

result <- data.frame(t=NA, iz=NA, icomb=NA, DGDq2=NA, dDGDq2=NA, chi=NA, p=NA, includesys=NA, fitcorr=NA)

for(time in c(48, 52, 56)) {
    file <- paste0("/home/gross/Documents/heavymesons/data/cB211ab.07.64/cov_diag/t", time, "/outputDGammaDq2/DGammaDq2.bin")
    if (covfull) file <- paste0("/home/gross/Documents/heavymesons/data/cB211ab.07.64/th6_t", time, "/outputDGammaDq2/DGammaDq2.bin")
    data <- read_in_DGDq2(filename=file, write=FALSE)
for (iset in c(0)) {
    for (iz in c(0, 1, 2, 3)) {
        for(icomb in c(0, 1, 2, 3, 4)) {
#~             print(data[[2+iset]]$epsilons)
# maybe replace this with a call to indexDG?
            namesel <- paste0("idg", ((iz)*5 + icomb)*data[[2+iset]]$metadata$neps + (0:(data[[2+iset]]$metadata$neps-1)), "ieps", 0:(data[[2+iset]]$metadata$neps-1), "icomb", icomb, "iz", iz)
            title=paste("iset", iset, "iz", iz, "icomb", icomb, "tsnk", time)
            print(title)
            
            DGamma <- unlist(sapply(X=data[[2+iset]][namesel], FUN=getElement, name="DGDq2mean")[1, ], use.names=F)
            dDGamma_sys <- unlist(sapply(X=data[[2+iset]][namesel], FUN=getElement, name="sys")[1, ], use.names=F)
            dDGamma <- sqrt(unlist(sapply(X=data[[2+iset]][namesel], FUN=getElement, name="DGDq2sd")[1, ], use.names=F))
            DGammaboot <- sapply(X=data[[2+iset]][namesel], FUN=getElement, name="DGDq2boots")[1, ]
            
            bs <- array(unlist(DGammaboot, use.names=F), dim=c(data[[2+iset]][[namesel[1]]]$bootnumber[[1]], data[[2+iset]]$metadata$neps))
            DGamma <- apply(X=bs, MARGIN=2, FUN=mean)
            print(DGamma)
            
            if (fitcorr) fitresult <- try(bootstrap.nlsfit(x=data[[2+iset]]$epsilons, y=DGamma, bsamples=bs, fn=fnpar, par.guess=c(1, 1, 1), mask=seq(1, uplim), verbose=F, CovMatrix = NULL))
            if (!fitcorr) fitresult <- try(bootstrap.nlsfit(x=data[[2+iset]]$epsilons, y=DGamma, bsamples=bs, fn=fnpar, par.guess=c(1, 1, 1), mask=seq(1, uplim), verbose=F))
            if(!inherits(fitresult, "try-error")) {
            try(plot(fitresult, plot.range=c(0, max(fitresult$x[seq(1, uplim)])), main=paste(title, "only stat"), xlab="epsilon", ylab="dG/dq2"))
            result <- rbind(result, data.frame(t=time, iz=iz, icomb=icomb, DGDq2=fitresult$t0[1], dDGDq2=fitresult$se[1], chi=fitresult$chisqr / fitresult$dof, p=fitresult$Qval, includesys=F, fitcorr=fitcorr))
            } else {
                result <- rbind(result, data.frame(t=time, iz=iz, icomb=icomb, DGDq2=NA, dDGDq2=NA, chi=NA, p=NA, includesys=NA, fitcorr=NA))
            }
            
            if(dosys) {
                bootsys <- parametric.bootstrap(boot.R=data[[2+iset]][[namesel[1]]]$bootnumber[[1]], x=rep(0, data[[2+iset]]$metadata$neps), dx=dDGamma_sys, seed=12345678)
                if (!fitcorr) fitresult <- try(bootstrap.nlsfit(x=data[[2+iset]]$epsilons, y=DGamma, bsamples=bs+bootsys, fn=fnpar, par.guess=c(1, 1, 1), mask=seq(1, uplim), CovMatrix=NULL))
                if (fitcorr) fitresult <- try(bootstrap.nlsfit(x=data[[2+iset]]$epsilons, y=DGamma, bsamples=bs+bootsys, fn=fnpar, par.guess=c(1, 1, 1), mask=seq(1, uplim)))
                if(!inherits(fitresult, "try-error")) {
#~                 print(fitresult)
                try(plot(fitresult, plot.range=c(0, max(fitresult$x[seq(1, uplim)])), main=paste(title, "total error"), xlab="epsilon", ylab="dG/dq2"))
                result <- rbind(result, data.frame(t=time, iz=iz, icomb=icomb, DGDq2=fitresult$t0[1], dDGDq2=fitresult$se[1], chi=fitresult$chisqr / fitresult$dof, p=fitresult$Qval, includesys=T, fitcorr=fitcorr))
                } else {
                    result <- rbind(result, data.frame(t=time, iz=iz, icomb=icomb, DGDq2=NA, dDGDq2=NA, chi=NA, p=NA, includesys=NA, fitcorr=NA))
                }
            }
            
        }
    }
}
}
result <- result[-1, ]
result <- result[ ! (result$iz==0 & (result$icomb==3 | result$icomb==4)), ]
print(result)
print(warnings())

savename <- "result_extrapolation_sigma_covdiag"
if (covfull) savename <- "result_extrapolation_sigma_covfull"
if (fitcorr) savename <- paste0(savename, "fitcorr.csv")
if (!fitcorr) savename <- paste0(savename, "fituncr.csv")
write.table(result, savename, row.names=F, col.names=T)
}

if (F) {
covfull <- F
savename <- "result_extrapolation_sigma_covdiag"
if (covfull) savename <- "result_extrapolation_sigma_covfull"
if (fitcorr) savename <- paste0(savename, "fitcorr.csv")
if (!fitcorr) savename <- paste0(savename, "fituncr.csv")

result <- read.table(savename, header=TRUE)
print(result)

newinfo <- data.frame(time=NA, icomb=NA, sumz0z1z2=NA, z3=NA, diffsumz3=NA, reldiff=NA)

for (time in unique(result$t)) {
    for (icomb in unique(result$icomb)) {
        sumresult <- sum(result$DGDq2[result$t==time & result$icomb==icomb & result$iz != 3])
        z3 <- result$DGDq2[result$t==time & result$icomb==icomb & result$iz == 3]
        newline <- data.frame(time=time, icomb=icomb, sumz0z1z2=sumresult, z3=z3, diffsumz3 = sumresult - z3, reldiff= 2*(sumresult - z3)/(sumresult + z3))
        newinfo <- rbind(newinfo, newline)
    }
}
newinfo <- na.omit(newinfo)
print(newinfo)
}
