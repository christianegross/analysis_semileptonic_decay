#' Do a complex division cf1 / cf2, taking real and imaginary part into account
#' (a+bi) / (c+di) = ((a+bi)*(c-di)) / (c^2 + d^2) = ((ac-bd) + (bc-ad)i) / (c^2 + d^2)
complexdivisionvlocal <- function(cf1, cf2){
  stopifnot(inherits(cf1, 'cf_meta'))
  stopifnot(inherits(cf2, 'cf_meta'))
  cf <- cf1
  if(inherits(cf1, 'cf_orig') && inherits(cf2, 'cf_orig')) {
    stopifnot(all(dim(cf1$cf) == dim(cf2$cf)))
    stopifnot(cf1$Time == cf2$Time)
    
    
    if( has_icf(cf1) | has_icf(cf2) ){
      stopifnot(has_icf(cf1) & has_icf(cf2))
      divisor <- cf2$cf^2 + cf2$icf^2
      # print(head(divisor[, 1:5]))
      # print(head(cf1$cf[, 1:5] * cf2$cf[, 1:5]))
      # print(head(cf1$icf[, 1:5] * cf2$icf[, 1:5]))
      cf$cf <- (cf1$cf * cf2$cf - cf1$icf * cf2$icf) / divisor
      cf$icf <- (cf1$icf * cf2$cf - cf1$cf * cf2$icf) / divisor
      ## the following is a bit dangerous, however, for
      ## principal correlators this is the only way to
      ## build a ratio
      if(inherits(cf1, 'cf_boot') && inherits(cf2, 'cf_boot') &&
         all(dim(cf1$cf.tsboot$t) == dim(cf2$cf.tsboot$t)) &&
         cf1$seed == cf2$seed && cf1$boot.l == cf2$boot.l) {
        
        divisor <- cf2$cf.tsboot$t^2 + cf2$icf.tsboot$t^2
        divisor0 <- cf2$cf.tsboot$t0^2 + cf2$icf.tsboot$t0^2
        
        cf$cf.tsboot$t  <- (cf1$cf.tsboot$t * cf2$cf.tsboot$t - cf1$icf.tsboot$t * cf2$icf.tsboot$t) / divisor
        cf$cf.tsboot$t0  <- (cf1$cf.tsboot$t0 * cf2$cf.tsboot$t0 - cf1$icf.tsboot$t0 * cf2$icf.tsboot$t0) / divisor0
        
        cf$icf.tsboot$t  <- (cf1$cf.tsboot$t * cf2$icf.tsboot$t - cf1$cf.tsboot$t * cf2$icf.tsboot$t) / divisor
        cf$icf.tsboot$t0  <- (cf1$cf.tsboot$t0 * cf2$icf.tsboot$t0 - cf1$cf.tsboot$t0 * cf2$icf.tsboot$t0) / divisor0
        
        cf$tsboot.se <- apply(cf$cf.tsboot$t, MARGIN = 2L, FUN = cf$error_fn)
        cf$cf0 <- cf$cf.tsboot$t0
        cf$itsboot.se <- apply(cf$icf.tsboot$t, MARGIN = 2L, FUN = cf$error_fn)
        cf$icf0 <- cf$icf.tsboot$t0
      }
    } else {
      cf <- cf1 / cf2
    }
  }
  else cf <- invalidate.samples.cf(cf)
  return (cf)
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
    if(inherits(cf1, 'cf_orig') && inherits(cf2, 'cf_orig')) {
      stopifnot(all(dim(cf1$cf) == dim(cf2$cf)))
      stopifnot(cf1$Time == cf2$Time)
      
      cf$cf <- cf1$cf / cf2$cf
      
      if( has_icf(cf1) | has_icf(cf2) ){
        stopifnot(has_icf(cf1) & has_icf(cf2))
        cf$icf <- cf1$icf / cf2$cf
      }
    }
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
    else cf <- invalidate.samples.cf(cf)
    # print(c(inherits(cf1, 'cf_boot'), inherits(cf2, 'cf_boot'),
    #   all(dim(cf1$cf.tsboot$t) == dim(cf2$cf.tsboot$t)),
    #   cf1$seed == cf2$seed, cf1$boot.l == cf2$boot.l))
    return (cf)
  }