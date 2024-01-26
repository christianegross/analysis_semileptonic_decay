
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


multiplycfbootstrap <- function(cf, central, bootstraps){
  cfnew <- cf
  cfnew$cf <- cf$cf * central
  if( has_icf(cf)){
    cfnew$icf <- cf$icf * central
  }
  if(inherits(cf, 'cf_boot')){
    # print(dim(cf$cf.tsboot$t))
    # print(dim(bootstraps))
    #' bootstraps may be an array, so use as.vector
    cfnew$cf.tsboot$t  <- cf$cf.tsboot$t * as.vector(bootstraps)
    # print(dim(cfnew$cf.tsboot$t))
    cfnew$cf.tsboot$t0  <- cf$cf.tsboot$t0 * central
    cfnew$tsboot.se <- apply(cfnew$cf.tsboot$t, MARGIN = 2L, FUN = cfnew$error_fn)
    cfnew$cf0 <- cfnew$cf.tsboot$t0
    
    if( has_icf(cf)){
      cfnew$icf.tsboot$t  <- cf$icf.tsboot$t * as.vector(bootstraps)
      cfnew$icf.tsboot$t0  <- cf$icf.tsboot$t0 * central
      
      cfnew$itsboot.se <- apply(cfnew$icf.tsboot$t, MARGIN = 2L, FUN = cfnew$error_fn)
      cfnew$icf0 <- cfnew$icf.tsboot$t0
    }
  }
  return(cfnew)
}
  
multiplycfscalar <- function(cf, scalar){
#~     print(scalar)
    cfnew <- cf
    cfnew$cf <- cf$cf * scalar
    if( has_icf(cf)){
      cfnew$icf <- cf$icf * scalar
    }
    if(inherits(cf, 'cf_boot')){
      # print(dim(cf$cf.tsboot$t))
      # print(dim(bootstraps))
      #' bootstraps may be an array, so use as.vector
      cfnew$cf.tsboot$t  <- cf$cf.tsboot$t * scalar
      # print(dim(cfnew$cf.tsboot$t))
      cfnew$cf.tsboot$t0  <- cf$cf.tsboot$t0 * scalar
      cfnew$tsboot.se <- apply(cfnew$cf.tsboot$t, MARGIN = 2L, FUN = cfnew$error_fn)
      cfnew$cf0 <- cfnew$cf.tsboot$t0
      
      if( has_icf(cf)){
        cfnew$icf.tsboot$t  <- cf$icf.tsboot$t * scalar
        cfnew$icf.tsboot$t0  <- cf$icf.tsboot$t0 * scalar
        
        cfnew$itsboot.se <- apply(cfnew$icf.tsboot$t, MARGIN = 2L, FUN = cfnew$error_fn)
        cfnew$icf0 <- cfnew$icf.tsboot$t0
      }
    }
    return(cfnew)
}

dividebyreal.cf <- function(cf1, cf2) {
    stopifnot(inherits(cf1, 'cf_meta'))
    stopifnot(inherits(cf2, 'cf_meta'))
    cf <- cf1
#~     print("assignment")
#~     print(names(cf))
    if(inherits(cf1, 'cf_orig') && inherits(cf2, 'cf_orig')) {
      stopifnot(all(dim(cf1$cf) == dim(cf2$cf)))
      stopifnot(cf1$Time == cf2$Time)
      
      cf$cf <- cf1$cf / cf2$cf
      
      if( has_icf(cf1) | has_icf(cf2) ){
        stopifnot(has_icf(cf1) & has_icf(cf2))
        cf$icf <- cf1$icf / cf2$cf
      }
    }
#~     print("before boot addition")
#~     print(names(cf))
    ## the following is a bit dangerous, however, for
    ## principal correlators this is the only way to
    ## build a ratio
    if(inherits(cf1, 'cf_boot') && inherits(cf2, 'cf_boot') &&
       all(dim(cf1$cf.tsboot$t) == dim(cf2$cf.tsboot$t)) &&
       cf1$seed == cf2$seed && cf1$boot.l == cf2$boot.l) {
      
      cf$cf.tsboot$t  <- cf1$cf.tsboot$t / cf2$cf.tsboot$t
      cf$cf.tsboot$t0 <- cf1$cf.tsboot$t0 / cf2$cf.tsboot$t0
      cf$tsboot.se <- apply(cf$cf.tsboot$t, MARGIN = 2L, FUN = cf$error_fn)
      cf$cf0 <- cf$cf.tsboot$t0
      if( has_icf(cf1) | has_icf(cf2) ){
        stopifnot(has_icf(cf1) & has_icf(cf2))
        cf$icf.tsboot$t  <- cf1$icf.tsboot$t / cf2$cf.tsboot$t
        cf$icf.tsboot$t0 <- cf1$icf.tsboot$t0 / cf2$cf.tsboot$t0
        cf$itsboot.se <- apply(cf$icf.tsboot$t, MARGIN = 2L, FUN = cf$error_fn)
        cf$icf0 <- cf$icf.tsboot$t0
      }
    }
    else {
#~         print("invalidate samples")
        cf <- invalidate.samples.cf(cf)
    }
    
#~     print("end")
#~     print(names(cf))
    # print(c(inherits(cf1, 'cf_boot'), inherits(cf2, 'cf_boot'),
    #   all(dim(cf1$cf.tsboot$t) == dim(cf2$cf.tsboot$t)),
    #   cf1$seed == cf2$seed, cf1$boot.l == cf2$boot.l))
    return (cf)
}

#~ void jack_store(FILE *ofp,jack_t *thisjack)
#~ {
#~    fwrite(&(thisjack->boot),sizeof(int),1,ofp);
#~    fwrite(&(thisjack->n),sizeof(int),1,ofp);
#~    fwrite(thisjack->f,sizeof(double),thisjack->n+4,ofp);
# bootstrapsamples, mean, meanwobias, sd, bias
#~ }

## some functions to tore the R-objects analogously to Nazario's code


## for the effective masses, we only store the result, not the masses for the timeslices
store_bin_cf_effmass <- function(to.write, cfeffmass, endian){
    if (cfeffmass$cf$resampling_method == "bootstrap") writeBin(object=as.integer(1), con=to.write, endian=endian)
    if (cfeffmass$cf$resampling_method == "jackknife") writeBin(object=as.integer(0), con=to.write, endian=endian)
    writeBin(object=as.integer(cfeffmass$boot.R), con=to.write, endian=endian)
    writeBin(object=as.double(cfeffmass$effmassfit$t[, 1]), con=to.write, endian=endian)
    writeBin(object=as.double(mean(cfeffmass$effmassfit$t[, 1])), con=to.write, endian=endian) # mean=mean(bootstrap samples)
    writeBin(object=as.double(cfeffmass$effmassfit$t0[1]), con=to.write, endian=endian) # meanwobias=effmass from mean values of effmass(t)
    writeBin(object=as.double(cfeffmass$effmassfit$se[1]), con=to.write, endian=endian) # sd=sd from cf object, sd of bootstrap samples
    writeBin(object=as.double(mean(cfeffmass$effmassfit$t[, 1]) - cfeffmass$effmassfit$t0[1]), con=to.write, endian=endian) # bias= mean(boostrapsamples) - effmass from mean values of effmass(t)    
}

## for the more general cf-objects, in this case the y-correlators, we store them timeslice by timeslice
## we can select whether to store the imaginary or real part of the correlator, which can be useful if we want to calculate something that has benn multiplid by i
store_bin_cf <- function(to.write, cf, endian, time, imre="reel"){
    if (cf$resampling_method == "bootstrap") writeBin(object=as.integer(1), con=to.write, endian=endian)
    if (cf$resampling_method == "jackknife") writeBin(object=as.integer(0), con=to.write, endian=endian)
    writeBin(object=as.integer(cf$boot.R), con=to.write, endian=endian)
    if(imre == "reel") {
        writeBin(object=as.double(cf$cf.tsboot$t[, time]), con=to.write, endian=endian)
#~         print(length(cf$cf.tsboot$t[, time]))
#~         print(cf$cf.tsboot$t[, time])
#~         print(cf$boot.R)
#~ print(cf$cf0[time])
        writeBin(object=as.double(mean(cf$cf.tsboot$t[, time])), con=to.write, endian=endian) # mean=mean(bootstrap samples)
        writeBin(object=as.double(cf$cf0[time]), con=to.write, endian=endian) # meanwobias=effmass from mean values of effmass(t)
        writeBin(object=as.double(cf$tsboot.se[time]), con=to.write, endian=endian) # sd=sd from cf object, sd of bootstrap samples
        writeBin(object=as.double(mean(cf$cf.tsboot$t[, time]) - cf$cf0[time]), con=to.write, endian=endian) # bias= mean(boostrapsamples) - effmass from mean values of effmass(t) 
    }   
    if(imre == "imag") {
#~ print(cf$cf0[time])
        writeBin(object=as.double(cf$icf.tsboot$t[, time]), con=to.write, endian=endian)
        writeBin(object=as.double(mean(cf$icf.tsboot$t[, time])), con=to.write, endian=endian) # mean=mean(bootstrap samples)
        writeBin(object=as.double(cf$icf0[time]), con=to.write, endian=endian) # meanwobias=effmass from mean values of effmass(t)
        writeBin(object=as.double(cf$itsboot.se[time]), con=to.write, endian=endian) # sd=sd from cf object, sd of bootstrap samples
        writeBin(object=as.double(mean(cf$icf.tsboot$t[, time]) - cf$icf0[time]), con=to.write, endian=endian) # bias= mean(boostrapsamples) - effmass from mean values of effmass(t) 
    }   
}


# we can also store an arbitrary vector
# vector has to contain mean as a first element and the bootstrapsamples after that
store_bin_vector <- function(to.write, vector, endian, resampling_method){
    if (resampling_method == "bootstrap") writeBin(object=as.integer(1), con=to.write, endian=endian)
    if (resampling_method == "jackknife") writeBin(object=as.integer(0), con=to.write, endian=endian)
    len <- length(vector)
    writeBin(object=as.integer(len-1), con=to.write, endian=endian)
    writeBin(object=as.double(vector[2:len]), con=to.write, endian=endian)
    writeBin(object=as.double(mean(vector[2:len])), con=to.write, endian=endian) # mean=mean of bootstrap
    writeBin(object=as.double(vector[1]), con=to.write, endian=endian) # meanwobias=mean from original data
    writeBin(object=as.double(sd(vector[2:len])), con=to.write, endian=endian) # sd=sd  of bootstrap samples
    writeBin(object=as.double(mean(vector[2:len]) - vector[1]), con=to.write, endian=endian) # bias= mean(boostrapsamples) - mean o original data    
}

store_bin_zero <- function(to.write, endian, resampling_method, boot.R) {
    if (resampling_method == "bootstrap") writeBin(object=as.integer(1), con=to.write, endian=endian)
    if (resampling_method == "jackknife") writeBin(object=as.integer(0), con=to.write, endian=endian)
    writeBin(object=as.integer(boot.R), con=to.write, endian=endian)
    writeBin(object=as.double(rep(0, (boot.R+4))), con=to.write, endian=endian) # fill everything with zeros
}
