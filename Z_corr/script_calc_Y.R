
library("hadron")
source("morefunctions.R")

args <- commandArgs(trailingOnly = TRUE)
inputfile <- as.character(args[1])
opt <- read.table(inputfile, sep=",", header=TRUE)
opt <- opt[1, ]



# TODO: include optparse for command line parameters or make inputfile
tsink <- opt$tsink
tj <- opt$tj
bootl <- opt$bs
bootsamples <- opt$nboot
maskfilelist <- 1:198
bootstraptype <- opt$statistics

if (bootstraptype != "jackknife" && bootstraptype != "bootstrap") stop("bootstraptype is invalid")

t1mH <- opt$t1heavy
t2mH <- opt$t2heavy

t1mL <- opt$t1light
t2mL <- opt$t2light

t1mL0 <- opt$t1light0
t2mL0 <- opt$t2light

theta <- opt$momentumindex




# change ZA, ZV with regards to previuos script
ZV <- opt$ZV
#~ dZV <- 0.00024
dZV <- 0
ZA <- opt$ZA
#~ dZA <- 0.000024
dZA <- 0
L <- opt$L
NT <- opt$T
afm <- opt$afm
mpigev <- opt$mpigev

thetain <- c(0, 0, opt$theta)
w <- sqrt(sum(thetain^2)) * pi / L
what <- 1

seed <- opt$seed



#~ Read in data for vector currents.


#~ define function for arbitrary $\mu, \nu$.
#~ <!-- In Antonios note, he only mentions AA and VV, so take out mixed parts. -->
#~ AV and VA are only needed for $C_{12}$ and $C_{21}$, not for the rest, so build a switch for them.
#~ The correlators are weighted with renormalization factors.
#~ reformulate Antonios note:

#~ $Y^2(t)=\frac{M_{D_s}}{\exp(M_{D_s}\cdot dt)}\frac{\hat{w}_0^2}{C_2(t)}\left(Z_V^2C_4^{V_0V_0}(t) + Z_A^2C_4^{A_0A_0}(t)\right)$
##four point function is not symmetric, do not symmetrise, also take care of division
# sym: for drawing bootstrapsamples, argument for tsboot. See documentation of hadron::bootstrap.cf
## seed is not fixed 
## to save time: allow to calculate only same indices or only mixed indices
getcmunu <- function(filellist, Time=128, mu=0, nu=0, boot.R=100, 
                    print=FALSE, boot.l=10, ZA, ZV, dZA, dZV, 
                    theta = 0, factor = 1, 
                    sim="fixed", bootstraptype="bootstrap", 
                    same = T, mixed = T, seed=123456){
aa <- NA
vv <- NA
av <- NA
va <- NA
if (same) {
aa <- readnissatextcf(filelist, smear_combs_to_read=c(paste("mes_contr_H_C_Dth", theta, "_A_C_P_H_H_S_H", sep="")), Time=Time, 
        corrtype="general", symmetrise=FALSE, sym.vec=c(1), 
        combs_to_read=data.frame(op1_idx="C_H", op2_idx=paste("Dth", theta, "_A", nu, "_C_P_H_H_S_H", sep=""), spin_comb=paste("A", mu, "P5", sep="")))
if (bootstraptype == "bootstrap") aa <- bootstrap.cf(aa, boot.R=boot.R, boot.l=boot.l, sim=sim, seed=seed)
if (bootstraptype == "jackknife") aa <- jackknife.cf(aa, boot.l=boot.l)

vv <- readnissatextcf(filelist, smear_combs_to_read=c(paste("mes_contr_H_C_Dth", theta, "_V_C_P_H_H_S_H", sep="")), Time=Time, 
        corrtype="general", symmetrise=FALSE, sym.vec=c(1), 
        combs_to_read=data.frame(op1_idx="C_H", op2_idx=paste("Dth", theta, "_V", nu, "_C_P_H_H_S_H", sep=""), spin_comb=paste("V", mu, "P5", sep="")))
if (bootstraptype == "bootstrap") vv <- bootstrap.cf(vv, boot.R=boot.R, boot.l=boot.l, sim=sim, seed=seed)
if (bootstraptype == "jackknife") vv <- jackknife.cf(vv, boot.l=boot.l)
}

# print(dim(bootstrapZA))
# print(dim(aa$cf.tsboot$t))

  
if (mixed) {
av <- readnissatextcf(filelist, smear_combs_to_read=c(paste("mes_contr_H_C_Dth", theta, "_V_C_P_H_H_S_H", sep="")), Time=Time, 
        corrtype="general", symmetrise=FALSE, sym.vec=c(1), 
        combs_to_read=data.frame(op1_idx="C_H", op2_idx=paste("Dth", theta, "_V", nu, "_C_P_H_H_S_H", sep=""), spin_comb=paste("A", mu, "P5", sep="")))
if (bootstraptype == "bootstrap") av <- bootstrap.cf(av, boot.R=boot.R, boot.l=boot.l, sim=sim, seed=seed)
if (bootstraptype == "jackknife") av <- jackknife.cf(av, boot.l=boot.l)

va <- readnissatextcf(filelist, smear_combs_to_read=c(paste("mes_contr_H_C_Dth", theta, "_A_C_P_H_H_S_H", sep="")), Time=Time, 
        corrtype="general", symmetrise=FALSE, sym.vec=c(1), 
        combs_to_read=data.frame(op1_idx="C_H", op2_idx=paste("Dth", theta, "_A", nu, "_C_P_H_H_S_H", sep=""), spin_comb=paste("V", mu, "P5", sep="")))
if (bootstraptype == "bootstrap") va <- bootstrap.cf(va, boot.R=boot.R, boot.l=boot.l, sim=sim, seed=seed)
if (bootstraptype == "jackknife") va <- jackknife.cf(va, boot.l=boot.l)
}


bootstrapZA <- parametric.bootstrap(boot.R, c(ZA), c(dZA), seed=seed)
bootstrapZV <- parametric.bootstrap(boot.R, c(ZV), c(dZV), seed=seed)


if (same) {
aa <- multiplycfbootstrap(cf = aa, central = ZA^2 * factor, bootstraps = bootstrapZA^2 * factor)
vv <- multiplycfbootstrap(cf = vv, central = ZV^2 * factor, bootstraps = bootstrapZV^2 * factor)
}

if (mixed) {
av <- multiplycfbootstrap(cf = av, central = ZA*ZV * factor, bootstraps = bootstrapZA*bootstrapZV * factor)
va <- multiplycfbootstrap(cf = va, central = ZV*ZA * factor, bootstraps = bootstrapZV*bootstrapZA * factor)
}

return(list(aa=aa, av=av, va=va, vv=vv))
}



## TODO: find a way to automatically exclude directories which do not have all files

#~ Make a filelist
filelist <- getorderedfilelist(path="/home/gross/Documents/heavymesons/data/th0/out_th0/", 
        basename="b", ending="/mes_contr_H_C_S_H")
filelist2 <- getorderedfilelist(path="/home/gross/Documents/heavymesons/data/th0/out_th0/", 
        basename="a", ending="/mes_contr_H_C_S_H")
filelist <- c(filelist, filelist2)

filelist <- substr(filelist, 1, nchar(filelist[1])-17)
#~ filelist <- filelist[-1]
filelist <- filelist[maskfilelist]



#~ symmetrisation makes 65 timeslices out of 128: measure from 0..127, 0 stays the same, 1+127, 2+126, ..., 64 stays the same -> now from 0..64, 65 times. Do not symmetrise four-point function, symmetrisation does not make sense with more than one inserted time. But do symmetrise two-point-function to get better signal for mass, matrix element, unsymmetrise later to get values for all timeslices. Calculate smeared-smeared contribution to determine masses.

mesonheavy <- readnissatextcf(filelist, smear_combs_to_read=c("/mes_contr_H_C_H_H_S_H"), Time=opt$T,
                corrtype="general", combs_to_read=data.frame(op1_idx="C_H", op2_idx="H_H_S_H", spin_comb="P5P5"),
                symmetrise=TRUE, sym.vec=c(1))
                
if (bootstraptype == "bootstrap") mesonheavy <- bootstrap.cf(mesonheavy, boot.R=bootsamples, boot.l=bootl, sim="fixed", seed=seed)
if (bootstraptype == "jackknife") mesonheavy <- jackknife.cf(mesonheavy, boot.l=bootl)
mesonheavyalltimeslices <- unsymmetrise.cf(mesonheavy)


mesonheavyeffmass <- bootstrap.effectivemass(mesonheavy, type="solve")
mesonheavyeffmass <- fit.effectivemass(mesonheavyeffmass, t1=t1mH, t2=t2mH, useCov=FALSE)
summary(mesonheavyeffmass)
#~ stop("")

w <- w/mesonheavyeffmass$effmassfit$t0[1]
w2 <- w^2
w2
bootsamples <- mesonheavy$boot.R
print(mesonheavy$resampling_method)

 # Contraction of Dth0_H ^ \dag and H_H_S_H
mesonlight <- readnissatextcf(filelist, smear_combs_to_read=c(paste("mes_contr_H_Dth", theta, "_H_H_S_H", sep="")), Time=opt$T,
                corrtype="general", combs_to_read=data.frame(op1_idx=paste0("Dth", theta, "_H"), op2_idx="H_H_S_H", spin_comb="P5P5"),
                symmetrise=TRUE, sym.vec=c(1))
if (bootstraptype == "bootstrap") mesonlight <- bootstrap.cf(mesonlight, boot.R=bootsamples, boot.l=bootl, sim="fixed", seed=seed)
if (bootstraptype == "jackknife") mesonlight <- jackknife.cf(mesonlight, boot.l=bootl)
mesonlightalltimeslices <- unsymmetrise.cf(mesonlight)


mesonlighteffmass <- bootstrap.effectivemass(mesonlight, type="solve")
mesonlighteffmass <- fit.effectivemass(mesonlighteffmass, t1=t1mL, t2=t2mL, useCov=FALSE)
summary(mesonlighteffmass)

 # Contraction of Dth0_H ^ \dag and H_H_S_H
mesonlight0 <- readnissatextcf(filelist, smear_combs_to_read=c(paste("mes_contr_H_Dth", theta, "_H_H_S_H", sep="")), Time=opt$T,
                corrtype="general", combs_to_read=data.frame(op1_idx=paste0("Dth", theta, "_H"), op2_idx="H_H_S_H", spin_comb="P5P5"),
                symmetrise=TRUE, sym.vec=c(1))
if (bootstraptype == "bootstrap") mesonlight0 <- bootstrap.cf(mesonlight0, boot.R=bootsamples, boot.l=bootl, sim="fixed", seed=seed)
if (bootstraptype == "jackknife") mesonlight0 <- jackknife.cf(mesonlight0, boot.l=bootl)
mesonlight0alltimeslices <- unsymmetrise.cf(mesonlight0)


mesonlight0effmass <- bootstrap.effectivemass(mesonlight0, type="solve")
mesonlight0effmass <- fit.effectivemass(mesonlight0effmass, t1=t1mL0, t2=t2mL0, useCov=FALSE)
summary(mesonlight0effmass)



#~ How to get $w_0$? Here, it is 1. For now, do not take ZA, ZV into account. At a later time, this can easily be changed.
#~ Determine four-point functions:
#~ From symmetry reasons, we get that we do not need to calculate all components, take them out with the mixed/same arguments
#~ We also do not need all matrix elements of cmunu

c_00 <- getcmunu(filellist = filelist, Time=opt$T, mu=0, nu=0, boot.R=bootsamples, ZA=ZA, ZV=ZV, dZA=dZA, dZV=dZV, factor=1, boot.l=bootl, same=T, mixed=F, bootstraptype=bootstraptype, seed=seed)
c_33 <- getcmunu(filellist = filelist, Time=opt$T, mu=3, nu=3, boot.R=bootsamples, ZA=ZA, ZV=ZV, dZA=dZA, dZV=dZV, factor=1, boot.l=bootl, same=T, mixed=F, bootstraptype=bootstraptype, seed=seed)
c_03 <- getcmunu(filellist = filelist, Time=opt$T, mu=0, nu=3, boot.R=bootsamples, ZA=ZA, ZV=ZV, dZA=dZA, dZV=dZV, factor=1, boot.l=bootl, same=T, mixed=F, bootstraptype=bootstraptype, seed=seed)
c_30 <- getcmunu(filellist = filelist, Time=opt$T, mu=3, nu=0, boot.R=bootsamples, ZA=ZA, ZV=ZV, dZA=dZA, dZV=dZV, factor=1, boot.l=bootl, same=T, mixed=F, bootstraptype=bootstraptype, seed=seed)
c_11 <- getcmunu(filellist = filelist, Time=opt$T, mu=1, nu=1, boot.R=bootsamples, ZA=ZA, ZV=ZV, dZA=dZA, dZV=dZV, factor=1, boot.l=bootl, same=T, mixed=F, bootstraptype=bootstraptype, seed=seed)
c_12 <- getcmunu(filellist = filelist, Time=opt$T, mu=1, nu=2, boot.R=bootsamples, ZA=ZA, ZV=ZV, dZA=dZA, dZV=dZV, factor=1, boot.l=bootl, same=F, mixed=T, bootstraptype=bootstraptype, seed=seed)
c_21 <- getcmunu(filellist = filelist, Time=opt$T, mu=2, nu=1, boot.R=bootsamples, ZA=ZA, ZV=ZV, dZA=dZA, dZV=dZV, factor=1, boot.l=bootl, same=F, mixed=T, bootstraptype=bootstraptype, seed=seed)
c_22 <- getcmunu(filellist = filelist, Time=opt$T, mu=2, nu=2, boot.R=bootsamples, ZA=ZA, ZV=ZV, dZA=dZA, dZV=dZV, factor=1, boot.l=bootl, same=T, mixed=F, bootstraptype=bootstraptype, seed=seed)





#~ <!-- is dt t_sink - t_j > 0  or t_j-t_sink < 0? -->
#~ <!-- here, $\hat{w}_0=1$. -->
#~ Determine prefactor for the H-tensors including bootstraps: $\frac{m_{D_s}}{\exp(m_{D_s} \cdot dt)}\cdot \frac{1}{C_2}$
#~ Take mass from bootstrapfit, to match the resampling procedure for $C_4$.
mds             <- c(mesonheavyeffmass$effmassfit$t0[1], mesonheavyeffmass$effmassfit$t[,1])
dt              <- tj-tsink
expmdsdt        <- exp(mds*dt)
mdsoverexpmdsdt <- mds / expmdsdt
expovermds      <- expmdsdt / mds
denom           <- multiplycfbootstrap(mesonheavyalltimeslices, expovermds[1], expovermds[2:length(mds)])


y_1 <- list(full=NA, vpar=NA, apar=NA, vperp=NA, vpar=NA)
y_2 <- list(full=NA, vpar=NA, apar=NA, vperp=NA, vpar=NA)
y_3 <- list(full=NA, vpar=NA, apar=NA, vperp=NA, vpar=NA)
y_4 <- list(full=NA, vpar=NA, apar=NA, vperp=NA, vpar=NA)
y_5 <- list(full=NA, vpar=NA, apar=NA, vperp=NA, vpar=NA)

## first fill the a, v, par, perp components according to the formula in the notes, then calculate the full as the sum of the components
## par: n=0
## perp: what=0

#~ y_1 <- c_11+c_22
y_1$vpar <- NA
y_1$apar <- NA
y_1$vperp <- c_11$vv+c_22$vv
y_1$aperp <- c_11$aa+c_22$aa
y_1$full <- y_1$vperp + y_1$aperp

#~ y_2 <- dividebyreal.cf(c_00, denom)
y_2$vpar <- dividebyreal.cf(c_00$vv, denom)
y_2$apar <- dividebyreal.cf(c_00$aa, denom)
y_2$aperp <- NA
y_2$vperp <- NA
y_2$full <- y_2$vpar + y_2$apar

#~ y_3 <- dividebyreal.cf(multiplycfscalar(c_33, -1), denom)
y_3$vpar <- dividebyreal.cf(multiplycfscalar(c_33$vv, -1*what*what), denom)
y_3$apar <- dividebyreal.cf(multiplycfscalar(c_33$aa, -1**what*what), denom)
y_3$aperp <- NA
y_3$vperp <- NA
y_3$full <- y_3$vpar + y_3$apar


#~ y_4 <- multiplycfscalar(c_03+c_30, 0.5)
y_4$apar <- multiplycfscalar(c_03$aa+c_30$aa, 0.5)
y_4$vpar <- multiplycfscalar(c_03$vv+c_30$vv, 0.5)
y_4$aperp <- NA
y_4$vperp <- NA
y_4$full <- y_4$vpar + y_4$apar

#~ y_5 <- c_12-c_21
y_5$vpar <- NA
y_5$apar <- NA
y_5$vperp <- NA
y_5$aperp <- NA
y_5$full <- (c_12$va+c_12$av)-(c_21$av+c_21$va)


for (name in c("full", "vpar", "apar", "vperp", "aperp")) {
    if (length(y_1[[name]]) > 1) y_1[[name]] <- dividebyreal.cf(multiplycfscalar(y_1[[name]], 0.5), denom)
    if (length(y_2[[name]]) > 1) y_2[[name]] <- multiplycfscalar(y_2[[name]], what*what)
    if (length(y_4[[name]]) > 1) y_4[[name]] <- dividebyreal.cf(multiplycfscalar(y_4[[name]], -1*what), denom)
    if (length(y_5[[name]]) > 1) y_5[[name]] <- dividebyreal.cf(multiplycfscalar(y_5[[name]], 0.5*what), denom)
}



#~ for(element in list(y_1, y_2, y_3, y_4, y_5)) {
#~     lapply(element, summary)
#~ }

## store
to.write <- file("Y.bin", "wb")
endian <- .Platform$endian

writeBin(object=as.integer(1),     con=to.write, endian=endian)
writeBin(object=as.integer(L),     con=to.write, endian=endian)
writeBin(object=as.integer(NT),    con=to.write, endian=endian)
writeBin(object=as.integer(tsink), con=to.write, endian=endian)
writeBin(object=as.integer(tj),    con=to.write, endian=endian)


writeBin(object=as.double(afm),    con=to.write, endian=endian)
writeBin(object=as.double(w),      con=to.write, endian=endian)
writeBin(object=as.double(mpigev), con=to.write, endian=endian)

store_bin_cf_effmass(to.write, mesonheavyeffmass,  endian)
store_bin_cf_effmass(to.write, mesonlighteffmass,  endian)
store_bin_cf_effmass(to.write, mesonlight0effmass, endian)

pdf("outputYbin/testY.pdf", title="")

ylist <- list(y_1, y_2, y_3, y_4, y_5)
# some Correlators are multiplied by i, so we have to store the imaginary part
imre <- c("reel", "reel", "reel", "imag", "imag")

for (comb in c("full", "vpar", "apar", "vperp", "aperp")) {
    for (index in seq(1, 5)) {
    y <- ylist[[index]]
    try(plot(y[[comb]], log="y", main=paste("y", index, "comb", comb), xlab="t", ylab="Y", xlim = c(0, 64)))
    try(plotwitherror(x=seq(0, 127), y=y[[comb]]$cf0, rep=TRUE, col=2))
        for (time in seq(tj+1, 1)) {
            if (length(y[[comb]]) > 1) {
                if (time==tsink) print(paste("y", index, "comb", comb))
                store_bin_cf(to.write, y[[comb]], endian, time=time, imre=imre[index])
            }
            else store_bin_zero(to.write, endian, resampling_method=bootstraptype, boot.R=bootsamples)
        }
    }
}

#~ int indexY(Y_t *Y,int t,int iy,int icomb)
#~ {
#~    return t+(Y->tj+1)*(iy+5*icomb);
#~ }
#~ {
#~    int iset,iy,ny;
   
#~    fwrite(&(nsets),sizeof(int),1,ofp);

#~    for(iset=0;iset<nsets;++iset)
#~    {
#~       fwrite(&(Y[iset].L),sizeof(int),1,ofp);
#~       fwrite(&(Y[iset].T),sizeof(int),1,ofp);
#~       fwrite(&(Y[iset].tsink),sizeof(int),1,ofp);
#~       fwrite(&(Y[iset].tj),sizeof(int),1,ofp);

#~       fwrite(&(Y[iset].afm),sizeof(double),1,ofp);
#~       fwrite(&(Y[iset].w),sizeof(double),1,ofp);
#~       fwrite(&(Y[iset].mpigev),sizeof(double),1,ofp);

#~       jack_store(ofp,Y[iset].mH);
#~       jack_store(ofp,Y[iset].mL);
#~       jack_store(ofp,Y[iset].eL);

#~       ny=5*NCOMBS*(Y[iset].tj+1);
#~       for(iy=0;iy<ny;++iy)
#~ 	 jack_store(ofp,Y[iset].y+iy);
#~    }
#~ }



#~ summary(y_1)
#~ summary(y_2)
#~ summary(y_3)
#~ summary(y_4)
#~ summary(y_5)

## TODO: write function to filter imaginary/real part
## TODO: write function to store it in bin
#~ data_christiane <- data.frame(t=seq(0,127), 
#~ c4=c_00$cf0, dc4=c_00$tsboot.se, 
#~ c2=mesonheavyalltimeslices$cf0, dc2=mesonheavyalltimeslices$tsboot.se)

#~ data_christiane$y_2 = c(y_2$cf0[(tj+1):1], rep(NA, 127-tj)) 
#~ data_christiane$dy2=c(y_2$tsboot.se[(tj+1):1], rep(NA, 127-tj))
#~ data_christiane$y_1 <- c(y_1$cf0[(tj+1):1], rep(NA, 127-tj))
#~ data_christiane$dy1 <- c(y_1$tsboot.se[(tj+1):1], rep(NA, 127-tj))
#~ data_christiane$y_3 <- c(y_3$cf0[(tj+1):1], rep(NA, 127-tj))
#~ data_christiane$dy3 <- c(y_3$tsboot.se[(tj+1):1], rep(NA, 127-tj))
#~ data_christiane$y_4 <- c(y_4$icf0[(tj+1):1], rep(NA, 127-tj))
#~ data_christiane$dy4 <- c(y_4$itsboot.se[(tj+1):1], rep(NA, 127-tj))
#~ data_christiane$y_5 <- c(y_5$icf0[(tj+1):1], rep(NA, 127-tj))
#~ data_christiane$dy5 <- c(y_5$itsboot.se[(tj+1):1], rep(NA, 127-tj))

#~ data_christiane$zzero <- data_christiane$y_2 + data_christiane$y_3 - 2 * data_christiane$y_4 
#~ data_christiane$dzzero <- sqrt(data_christiane$dy2^2 + data_christiane$dy3^2 + 4 * data_christiane$dy4^2)
#~ data_christiane$zone <- - 4 * data_christiane$y_1 + 2 * data_christiane$y_3 - 2 * data_christiane$y_4 
#~ data_christiane$dzone <- sqrt(16 * data_christiane$dy1^2 + 4 * data_christiane$dy3^2 + 4 * data_christiane$dy4^2)
#~ data_christiane$ztwo <- data_christiane$y_3 - 2 * data_christiane$y_1 
#~ data_christiane$dztwo <- sqrt(data_christiane$dy3^2 + 4 * data_christiane$dy1^2)
#~ write.table(data_christiane, "correlators_christiane_th0_tj48.csv", row.names=F, col.names=T)

#~ Summary of all Ys, with factors, im, real used
#~ ```{r, echo=FALSE, eval=TRUE}
#~ summary <- data.frame(
#~   corr = c("Y1", "Y2", "Y3", "Y4", "Y5"),
#~   indices = c("11+22", "00", "-33", "-1/2*(03+30)", "1/2*(12-21)"),
#~   realpart =  c("RE", "RE", "RE", "IM", "IM"), 
#~   AAplusVV = c(T, T, T, T, F),
#~   AVplusVA = c(F, F, F, F, T),
#~   additionalfactor = c(0.5, 1, 1, 1, 1)
#~ )
#~ knitr::kable(summary)
#~ ```

#~ Antonio wrote that the additional factor $\frac{1}{2}$ in $Y^1$ comes from the fact that the contribution is equally split between two unit vectors.  This is not the case for the other $Y$.
#~ The factor $\frac{1}{2}$ in $Y^4$ and $Y^5$ comes from the paper.
