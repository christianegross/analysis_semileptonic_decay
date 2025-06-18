library("hadron")
library("christianesfunctions")
# devtools::load_all("~/code/christianesfunctions")
args <- commandArgs(trailingOnly = TRUE)
# 1: input file
# 2: output folder

if(file.exists(sprintf("%s/Y.log", args[2]))) stop("Y output file already exists")
sink(sprintf("%s/Y.log", args[2]))


## read in analysis parameters
parameters <- read.table(args[1], fill = TRUE, col.names = paste0("V", 1:4))


## determine mass heavy
Hmass <- readcffrompreprocessed(realfile = paste(sep = "/", parameters$V2[parameters$V1 == "path"], parameters$V2[parameters$V1 == "Hfile"]), symmetrise = T)


if (parameters$V2[parameters$V1 == "statistics"] == "BOOT") Hmass <- bootstrap.cf(Hmass, boot.R = as.numeric(parameters$V2[parameters$V1 == "nboot"]), boot.l = as.numeric(parameters$V2[parameters$V1 == "bs"]))
if (parameters$V2[parameters$V1 == "statistics"] == "JACK") Hmass <- jackknife.cf(Hmass, boot.l = as.numeric(parameters$V2[parameters$V1 == "bs"]))

emass <- bootstrap.effectivemass(Hmass, type = "solve")

emass <- fit.effectivemass(emass, t1 = as.numeric(parameters$V2[parameters$V1 == "tfith"]), t2 = as.numeric(parameters$V3[parameters$V1 == "tfith"]))

summary(emass, verbose=F)
print(data.frame(t= seq_along(emass$t0), m = emass$t0, dm = emass$se))




## determine mass light
Lmass <- readcffrompreprocessed(realfile = paste(sep = "/", parameters$V2[parameters$V1 == "path"], parameters$V2[parameters$V1 == "Lfile"]), symmetrise = T)


if (parameters$V2[parameters$V1 == "statistics"] == "BOOT") Lmass <- bootstrap.cf(Lmass, boot.R = as.numeric(parameters$V2[parameters$V1 == "nboot"]), boot.l = as.numeric(parameters$V2[parameters$V1 == "bs"]))
if (parameters$V2[parameters$V1 == "statistics"] == "JACK") Lmass <- jackknife.cf(Lmass, boot.l = as.numeric(parameters$V2[parameters$V1 == "bs"]))

Lemass <- bootstrap.effectivemass(Lmass, type = "solve")

Lemass <- fit.effectivemass(Lemass, t1 = as.numeric(parameters$V2[parameters$V1 == "tfitl"]), t2 = as.numeric(parameters$V3[parameters$V1 == "tfitl"]))

## determine mass heavy
L0mass <- readcffrompreprocessed(realfile = paste(sep = "/", parameters$V2[parameters$V1 == "path"], parameters$V2[parameters$V1 == "Hfile"]), symmetrise = T)


if (parameters$V2[parameters$V1 == "statistics"] == "BOOT") L0mass <- bootstrap.cf(L0mass, boot.R = as.numeric(parameters$V2[parameters$V1 == "nboot"]), boot.l = as.numeric(parameters$V2[parameters$V1 == "bs"]))
if (parameters$V2[parameters$V1 == "statistics"] == "JACK") L0mass <- jackknife.cf(L0mass, boot.l = as.numeric(parameters$V2[parameters$V1 == "bs"]))

L0emass <- bootstrap.effectivemass(L0mass, type = "solve")

L0emass <- fit.effectivemass(L0emass, t1 = as.numeric(parameters$V2[parameters$V1 == "tfitl0"]), t2 = as.numeric(parameters$V3[parameters$V1 == "tfitl0"]))

## todo: determine mass and effectivemass also on original data
## or determine Y-correlators only on bootstrap samples

## read in correlators
correlators <- list()


for (firstindex in c("A", "V")) {
    for (secondindex in c("A", "V")) {
        for (i in 0:3) {
            for (j in 0:3) {
                correlators[[paste0(firstindex, i, secondindex, j)]] <- readcffrompreprocessed(realfile = paste(sep = "/", parameters$V2[parameters$V1 == "path"], paste0(firstindex, i, secondindex, j, ".dat")), symmetrise = F)
            }
        }
    }
}
if (parameters$V2[parameters$V1 == "statistics"] == "BOOT") correlators <- lapply(correlators, bootstrap.cf, boot.R = as.numeric(parameters$V2[parameters$V1 == "nboot"]), boot.l = as.numeric(parameters$V2[parameters$V1 == "bs"]))
if (parameters$V2[parameters$V1 == "statistics"] == "JACK") correlators <- lapply(correlators, jackknife.cf, boot.l = as.numeric(parameters$V2[parameters$V1 == "bs"]))

#~ invisible(lapply(correlators, summary))


## determine Y

## calculate w, what
characteristics <- list()


characteristics$what <- c(1, as.numeric(parameters[parameters$V1 == "theta", 2:4]))
characteristics$w <- sqrt(sum(characteristics$what[2:4]^2))
characteristics$what[2:4] <- characteristics$what[2:4] / characteristics$w
characteristics$w <- characteristics$w * pi / as.numeric(parameters$V2[parameters$V1 == "L"])
characteristics$w <- characteristics$w / emass$effmassfit$t0[1]
characteristics$n1 <- c(0, as.numeric(parameters[parameters$V1 == "n1", 2:4]))
characteristics$n2 <- c(0, as.numeric(parameters[parameters$V1 == "n2", 2:4]))

characteristics$za <- as.numeric(parameters$V2[parameters$V1 == "za"])
characteristics$zv <- as.numeric(parameters$V2[parameters$V1 == "zv"])
tj <- as.numeric(parameters$V2[parameters$V1 == "tj"])
tsink <- as.numeric(parameters$V2[parameters$V1 == "tsink"])

print(characteristics)
print(Hmass)


tmp <- unsymmetrise.cf(Hmass)
ret <- tmp
ret$cf.tsboot$t <- tmp$cf.tsboot$t * exp((tj - tsink) * emass$effmassfit$t[, 1]) / emass$effmassfit$t[, 1]
ret$cf.tsboot$t0 <- tmp$cf.tsboot$t0 * exp((tj - tsink) * emass$effmassfit$t0[1]) / emass$effmassfit$t0[1]
ret <- cf_boot(ret, boot.R=ret$boot.R, boot.l=ret$boot.l, seed=ret$seed, sim=ret$sim, endcorr=ret$endcorr, cf.tsboot=ret$cf.tsboot, resampling_method=ret$resampling_method)
# ret$cf0
#~ boot.R, boot.l, seed, sim, endcorr, cf.tsboot, icf.tsboot = NULL, resampling_method
#~ ret$tsboot.se <- apply(X = ret$cf.tsboot$t, MARGIN = 2, FUN = sd)
#~ ret$cf0 <- ret$cf.tsboot$t0
#~ exp((tj - tsink) * emass$effmassfit$t0[1]) / emass$effmassfit$t0[1]
#~ emass$effmassfit$t0[1]

#   After we calculate the y, we have to adjust the time indexing. We do this slice-by-slice.

setcorrecttimeindex <- function(y, tj) {
    tmp <- y
    y$cf.tsboot$t[, (0:(tj)) + 1] <- tmp$cf.tsboot$t[, tj - (0:(tj)) + 1]
    y$cf.tsboot$t0[(0:(tj)) + 1] <- tmp$cf.tsboot$t0[tj - (0:(tj)) + 1]
    y$cf.tsboot$t[, (tj + 2):length(y$cf.tsboot$t[1, ])] <- 0
    y$cf.tsboot$t0[(tj + 2):length(y$cf.tsboot$t[1, ])] <- 0
    y <- cf_boot(y, boot.R=y$boot.R, boot.l=y$boot.l, seed=y$seed, sim=y$sim, endcorr=y$endcorr, cf.tsboot=y$cf.tsboot, resampling_method=y$resampling_method)
#~     y$tsboot.se <- apply(X = y$cf.tsboot$t, MARGIN = 2, FUN = sd)
#~     y$cf0 <- y$cf.tsboot$t0
    return(y)
}

yallcomb <- list()
for (comb in c("FULL", "VPAR", "APAR", "VPERP", "APERP")) {
    cat(paste("## Now analysing combination", comb), sep="\n")
    y1 <- zero.cf(nsamples = Hmass$conf.index, hasim = F, doboot = T, nrObs = 1, Time = ret$Time, nrStypes = 1, symmetrised = F, boot.R = ret$boot.R, boot.l = ret$boot.l, seed = ret$seed, sim = ret$sim, endcorr = ret$endcorr, resampling_method = ret$resampling_method)
    y2 <- zero.cf(nsamples = Hmass$conf.index, hasim = F, doboot = T, nrObs = 1, Time = ret$Time, nrStypes = 1, symmetrised = F, boot.R = ret$boot.R, boot.l = ret$boot.l, seed = ret$seed, sim = ret$sim, endcorr = ret$endcorr, resampling_method = ret$resampling_method)
    y3 <- zero.cf(nsamples = Hmass$conf.index, hasim = F, doboot = T, nrObs = 1, Time = ret$Time, nrStypes = 1, symmetrised = F, boot.R = ret$boot.R, boot.l = ret$boot.l, seed = ret$seed, sim = ret$sim, endcorr = ret$endcorr, resampling_method = ret$resampling_method)
    y4 <- zero.cf(nsamples = Hmass$conf.index, hasim = F, doboot = T, nrObs = 1, Time = ret$Time, nrStypes = 1, symmetrised = F, boot.R = ret$boot.R, boot.l = ret$boot.l, seed = ret$seed, sim = ret$sim, endcorr = ret$endcorr, resampling_method = ret$resampling_method)
    y5 <- zero.cf(nsamples = Hmass$conf.index, hasim = F, doboot = T, nrObs = 1, Time = ret$Time, nrStypes = 1, symmetrised = F, boot.R = ret$boot.R, boot.l = ret$boot.l, seed = ret$seed, sim = ret$sim, endcorr = ret$endcorr, resampling_method = ret$resampling_method)

    what <- characteristics$what
    n1 <- characteristics$n1
    n2 <- characteristics$n2
    za <- characteristics$za
    zv <- characteristics$zv


    if (comb == "APAR" || comb == "APERP") zv <- 0.0
    if (comb == "VPAR" || comb == "VPERP") za <- 0.0

    if (comb == "APAR" || comb == "VPAR") {
        n1 <- c(0.0, 0.0, 0.0, 0.0)
        n2 <- c(0.0, 0.0, 0.0, 0.0)
    }
    if (comb == "APERP" || comb == "VPERP") what <- c(0.0, 0.0, 0.0, 0.0)

    ## Y1
    cat("## Y1", sep="\n")
    for (i in 1:3) {
        for (j in 1:3) {
            tmp <- mul.cf(a = zv^2, cf = correlators[[paste0("V", i, "V", j)]])
            tmp <- tmp + mul.cf(a = za^2, cf = correlators[[paste0("A", i, "A", j)]])
            tmp <- mul.cf(cf = tmp, a = 0.5 * (n1[i + 1] * n1[j + 1] + n2[i + 1] * n2[j + 1]))
            y1 <- y1 + tmp
        }
    }

    y1 <- y1 / ret
    y1 <- setcorrecttimeindex(y = y1, tj = tj)


    ## Y2
    cat("## Y2")
    y2 <- mul.cf(a = zv^2, cf = correlators[[paste0("V", 0, "V", 0)]]) + mul.cf(a = za^2, cf = correlators[[paste0("A", 0, "A", 0)]])
    plot(y2, log = "y")
    plot(ret, log = "y")
    y2 <- y2 / ret
    y2 <- setcorrecttimeindex(y2, tj = tj)


    ## Y3
    cat("## Y3", sep="\n")

    for (i in 1:3) {
        for (j in 1:3) {
            tmp <- mul.cf(a = zv^2, cf = correlators[[paste0("V", i, "V", j)]]) + mul.cf(a = za^2, cf = correlators[[paste0("A", i, "A", j)]])
            tmp <- mul.cf(a = -what[i + 1] * what[j + 1], cf = tmp)
            y3 <- y3 + tmp
        }
    }
    y3 <- y3 / ret
    y3 <- setcorrecttimeindex(y3, tj = tj)

    # plot(y3, log = "y")



    ## Y4
    cat("## Y4", sep="\n")
    for (i in 1:3) {
        tmp <- mul.cf(a = zv^2, cf = correlators[[paste0("V", i, "V", 0)]]) + mul.cf(a = za^2, cf = correlators[[paste0("A", i, "A", 0)]])
        tmp2 <- mul.cf(a = zv^2, cf = correlators[[paste0("V", 0, "V", i)]]) + mul.cf(a = za^2, cf = correlators[[paste0("A", 0, "A", i)]])
        tmp <- tmp + tmp2
        tmp <- mul.cf(a = -0.5 * what[i + 1], cf = tmp)
        y4 <- y4 + tmp
    }
    y4 <- y4 / ret
    y4 <- setcorrecttimeindex(y4, tj = tj)


    ## Y5
    
    cat("## Y5", sep="\n")
    for (k in 1:3) {
        i <- k %% 3 + 1
        j <- (k + 1) %% 3 + 1
        tmpfactor <- 0.5 * what[k + 1] * za * zv
        tmp <- mul.cf(a = +tmpfactor, correlators[[paste0("V", i, "A", j)]] + correlators[[paste0("A", i, "V", j)]])
        tmp2 <- mul.cf(a = -tmpfactor, correlators[[paste0("V", j, "A", i)]] + correlators[[paste0("A", j, "V", i)]])

        y5 <- y5 + (tmp + tmp2)
    }
    y5 <- y5 / ret
    y5 <- setcorrecttimeindex(y5, tj = tj)

    yallcomb[[paste0("y", 1, comb)]] <- c(y1)
    yallcomb[[paste0("y", 2, comb)]] <- c(y2)
    yallcomb[[paste0("y", 3, comb)]] <- c(y3)
    yallcomb[[paste0("y", 4, comb)]] <- c(y4)
    yallcomb[[paste0("y", 5, comb)]] <- c(y5)
}

## TODO: store results, properly document functions in storebin


## store
cat("## Now storing results", sep="\n")
to.write <- file(sprintf("%s/Y.bin", args[2]), "wb")
endian <- .Platform$endian

writeBin(object=as.integer(1),                                       con=to.write, endian=endian)
writeBin(object=as.integer(parameters$V2[parameters$V1 == "L"]),     con=to.write, endian=endian)
writeBin(object=as.integer(parameters$V2[parameters$V1 == "T"]),     con=to.write, endian=endian)
writeBin(object=as.integer(characteristics$tsink),                   con=to.write, endian=endian)
writeBin(object=as.integer(characteristics$tj),                      con=to.write, endian=endian)


writeBin(object=as.double(parameters$V2[parameters$V1 == "afm"]),    con=to.write, endian=endian)
writeBin(object=as.double(characteristics$w),                        con=to.write, endian=endian)
writeBin(object=as.double(parameters$V2[parameters$V1 == "mpigev"]), con=to.write, endian=endian)

store_bin_cf_effmass(to.write, emass,   endian)
store_bin_cf_effmass(to.write, Lemass,  endian)
store_bin_cf_effmass(to.write, L0emass, endian)

for (comb in c("FULL", "VPAR", "APAR", "VPERP", "APERP")) {
    for(index in 1:5) {
        for(time in seq(1, tj+1)) {
            if(length(yallcomb[[paste0("y", index, comb)]]$cf0) > 1) {
            store_bin_cf(to.write, yallcomb[[paste0("y", index, comb)]], endian, time=time)
            } else {
                store_bin_zero(to.write, endian, resampling_method=ifelse((parameters$V2[parameters$V1 == "statistics"] == "BOOT"), "bootstrap", "jackknife"), boot.R=ret$boot.R)
            }
        }
    }
}
close(to.write)

cat("## Stored results", sep="\n")

sink(file=NULL)
sink(sprintf("%s/Y.dat", args[2]))
cat("# $1  mh", "# $2  dmh", "# $3  iset", "# $4  w", "# $5  iy", "# $6  t", "# $7  Y_FULL", "# $8  dY_FULL", "# $9  Y_VPAR", "# $10 dY_VPAR", "# $11 Y_APAR", "# $12 dY_APAR", "# $13 Y_VPERP", "# $14 dY_VPERP", "# $15 Y_APERP", "# $16 dY_APERP", sep="\n")
for(index in 1:5) {
    yres <- data.frame(mass=emass$effmassfit$t0[1], dmass=emass$effmassfit$se[1], iset=0, w=characteristics$w, iy=index, t=(0:(tj)))
    for (comb in c("FULL", "VPAR", "APAR", "VPERP", "APERP")) {
        yres[, paste0(comb, "mn")] <- yallcomb[[paste0("y", index, comb)]]$cf0[(0:(tj))+1]
        yres[, paste0(comb, "sd")] <- yallcomb[[paste0("y", index, comb)]]$tsboot.se[(0:(tj))+1]
    }
    
    for(i in 1:dim(yres)[1]){
        cat(sprintf("%20.6e %20.6e %4d %20.6e %d %4d %20.6e %20.6e %20.6e %20.6e %20.6e %20.6e %20.6e %20.6e %20.6e %20.6e ", 
        yres[i, 1], yres[i, 2], yres[i, 3], yres[i, 4], yres[i, 5], yres[i, 6], yres[i, 7], yres[i, 8], yres[i, 9], yres[i, 10], yres[i, 11], yres[i, 12], yres[i, 13], yres[i, 14], yres[i, 15], yres[i, 16]))
        cat("\n")
    }
    cat("\n\n")
}

sink(file=NULL)

# print("hi")

# for(index in 1:5) {
    # print(yallcomb[[paste0("y", index, "FULL")]]$cf0)
# }


warnings()
#~ stop("hello")
