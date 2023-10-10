library("hadron")
source("read_in_binary.R")

#~ data <- read_in_DGDq2(filename="/home/gross/Documents/heavymesons/data/cB211ab.07.64/t48_64/outputDGammaDq2/DGammaDq2.bin", write=FALSE)
#~ data <- read_in_DGDq2(filename="/home/gross/Documents/heavymesons/data/cB211ab.07.64/th6_t48/outputDGammaDq2_10sigma/DGammaDq2.bin", write=FALSE)


fnlin <- function(par, x, boot.R, ...) par[1] + par[2] * x
fncon <- function(par, x, boot.R, ...) par[1]

filenames <- paste0("/home/gross/Documents/heavymesons/data/cB211ab.07.64/th6_t", c(48, 52, 56), "/outputDGammaDq2_10sigma/DGammaDq2.bin")

pdf("extrapolate_epsilon.pdf", title="")
result <- data.frame(t=NA, iz=NA, icomb=NA, DGDq2=NA, dDGDq2=NA, chi=NA, p=NA)
for(time in c(48)) {
    file <- paste0("/home/gross/Documents/heavymesons/data/cB211ab.07.64/th6_t", time, "/outputDGammaDq2_10sigma/DGammaDq2.bin")
    data <- read_in_DGDq2(filename=file, write=FALSE)
for (iset in c(0)) {
    for (iz in seq(0, 3)) {
        for(icomb in seq(0, 4)) {
            print(data[[2+iset]]$epsilons)
            namesel <- paste0("idg", ((iz)*5 + icomb)*data[[2+iset]]$metadata$neps + (0:(data[[2+iset]]$metadata$neps-1)), "ieps", 0:(data[[2+iset]]$metadata$neps-1), "icomb", icomb, "iz", iz)
            title=paste("iset", iset, "iz", iz, "icomb", icomb, "tsnk", time)
            
            
    #~         DGamma <- c()
    #~         dDGamma <- c()
    #~         for (name in namesel) {
    #~             DGamma <- append(DGamma, data[[2]][[name]]$DGDq2mean[[1]][1])
    #~             dDGamma <- append(dDGamma, data[[2]][[name]]$DGDq2sd[[1]][1])
    #~         }
            DGamma <- unlist(sapply(X=data[[2+iset]][namesel], FUN=getElement, name="DGDq2mean")[1, ], use.names=F)
            dDGamma <- sqrt(unlist(sapply(X=data[[2+iset]][namesel], FUN=getElement, name="DGDq2sd")[1, ], use.names=F))
            DGammaboot <- sapply(X=data[[2+iset]][namesel], FUN=getElement, name="DGDq2boots")[1, ]
#~             print(names(data[[2+iset]][[namesel[1]]]))
#~             print(data[[2+iset]][[namesel[1]]]$bootnumber[[1]])
            bs <- array(unlist(DGammaboot, use.names=F), dim=c(data[[2+iset]][[namesel[1]]]$bootnumber[[1]], data[[2+iset]]$metadata$neps))
#~             print(bs)
            
            try(plotwitherror(x=data[[2+iset]]$epsilons, y=DGamma, dy=dDGamma, xlab="epsilon", ylab="DG/Dq2", main=title))
            
            try(fitresult <- bootstrap.nlsfit(x=data[[2+iset]]$epsilon, y=DGamma, bsamples=bs, fn=fnlin, par.guess=c(1, 1), mask=seq(1, 6)))
            try(print(fitresult))
            try(plot(fitresult, plot.range=c(0, max(fitresult$x[seq(1, 6)])), main=title, xlab="epsilon", ylab="dG/dq2"))
            
            fitresult <- try(bootstrap.nlsfit(x=data[[2+iset]]$epsilons, y=DGamma, bsamples=bs, fn=fncon, par.guess=c(1), mask=seq(1, 6)))
            if(!inherits(fitresult, "try-error")) {
            print(fitresult)
            try(plot(fitresult, plot.range=c(0, max(fitresult$x[seq(1, 6)])), main=title, xlab="epsilon", ylab="dG/dq2"))
            print(names(fitresult))
            result <- rbind(result, data.frame(t=time, iz=iz, icomb=icomb, DGDq2=fitresult$t0[1], dDGDq2=fitresult$se[1], chi=fitresult$chisqr / fitresult$dof, p=fitresult$Qval))
            } else {
                result <- rbind(result, data.frame(t=time, iz=iz, icomb=icomb, DGDq2=NA, dDGDq2=NA, chi=NA, p=NA))
            }
            
        }
    }
}
}
result <- result[-1, ]
print(result)
print(warnings())
#~ print(namesel)
#~ print(names(data[[2]][[namesel[1]]]))

#~ sapply(X=data[[2]][namesel], FUN=getElement, name="sys")
#~ res <- sapply(X=data[[2]][namesel], FUN=getElement, name="DGDq2mean")
#~ res
#~ str(res)
#~ as.vector(res)
#~ res[1, ]
#~ DGamma <- unlist(sapply(X=data[[2]][namesel], FUN=getElement, name="DGDq2mean")[1, ], use.names=F)
#~ DGammaboot <- sapply(X=data[[2]][namesel], FUN=getElement, name="DGDq2boots")[1, ]
#~ print(DGammaboot)
#~ print(unlist(DGammaboot, use.names=F))
#~ bs <- array(unlist(DGammaboot, use.names=F), dim=c(15, 10))
#~ dDGamma <- unlist(sapply(X=data[[2]][namesel], FUN=getElement, name="DGDq2sd")[1, ], use.names=F)
#~ names(DGamma) <- NULL
#~ unlist(DGamma, use.names=F)
#~ dDGamma
#~ try(plotwitherror(x=data[[2]]$epsilon, y=DGamma, dy=dDGamma, xlab="eps", ylab="DG/Dq2", main=paste("hi")))

fnlin <- function(par, x, boot.R, ...) par[1] + par[2] * x

#~ fitresult <- bootstrap.nlsfit(x=data[[2]]$epsilon, y=DGamma, bsamples=bs, fn=fnlin, par.guess=c(1, 1), mask=seq(1, 6))
#~ print(fitresult)
#~ plot(fitresult, plot.range=c(0, max(fitresult$x[seq(1, 6)])))

fncon <- function(par, x, boot.R, ...) par[1]

#~ fitresult <- bootstrap.nlsfit(x=data[[2]]$epsilon, y=DGamma, bsamples=bs, fn=fncon, par.guess=c(1), mask=seq(1, 6))
#~ print(fitresult)
#~ plot(fitresult)
