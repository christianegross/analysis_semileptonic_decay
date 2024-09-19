## implements the trapezoidal integration for arbitrary distances 
## between the supporting points
## Order of y, x: can put several bootstrapsamples with common x-coordinates with apply
## x: supporting point; y=f(x)=function value
trapezoidal <- function(yval, xval, continue=FALSE, higherlimit=0, replacelower=FALSE, lowerlimit=0, replaceindex=0) {
    stopifnot(!(continue && replacelower))
    stopifnot(length(xval)==length(yval))
    if(any(is.na(yval))) return(NA)
    # stopifnot(xval==sort(x=xval) || xval==rev(sort(x=xval)))
    sum <- 0
    
    if(!(continue || replacelower)) {
        for(i in seq(1, length(xval)-1)) {
            sum <- sum + 1/2*(xval[i+1]-xval[i])*(yval[i]+yval[i+1])
        }
    }
    if(continue) {
    ## calculate slope, add end values to xval, yval
        stopifnot(higherlimit > max(xval))
        len <- length(xval)
        slope <- (yval[len] - yval[len-1])/(xval[len] - xval[len-1])
        higher <- yval[len] + (higherlimit-xval[len]) * slope
        yval <- append(yval, higher)
        xval <- append(xval, higherlimit)
        for(i in seq(1, length(xval)-1)) {
            sum <- sum + 1/2*(xval[i+1]-xval[i])*(yval[i]+yval[i+1])
        }
    }
    if(replacelower) {
        stopifnot(replaceindex > 0 && replaceindex <= length(xval))
        xval <- xval[1:replaceindex]
        len <- length(xval)
        slope <- (yval[len] - yval[len-1])/(xval[len] - xval[len-1])
        higher <- yval[len] + (lowerlimit-xval[len]) * slope
        xval[replaceindex] <- lowerlimit
        yval[len] <- higher
        for(i in seq(1, length(xval)-1)) {
            sum <- sum + 1/2*(xval[i+1]-xval[i])*(yval[i]+yval[i+1])
        }
    }
    return(sum)
}

## simpson's rule for an odd number of points, points do not have to beequidistant
simpson_odd <- function(yval, xval) {
    stopifnot(length(xval)==length(yval))
    stopifnot(length(xval)%%2==1)
    stopifnot(xval==sort(x=xval) || xval==rev(sort(x=xval)))
    sum <- 0
    for(i in seq(1, length(xval)-2, 2)) {
        sum <- sum + 1/6*(xval[i+2]-xval[i])*
                (yval[i] * (-3*xval[i+1]+2*xval[i]+xval[i+2])/(xval[i]-xval[i+1])
                +yval[i+1] * (xval[i] - xval[i+2])^2 / ((xval[i+1]-xval[i])*(xval[i+2]-xval[i+1]))
                +yval[i+2] * (-3*xval[i+1]+2*xval[i+2]+xval[i])/(xval[i+2]-xval[i+1]))
    }
    return (sum)
}

## Generalisation of simpson's rule for odd and even number of points
## take trapezoidal rule for endpoints, average over both possibilities 
## (trapezoidal rule lower boundary or upper boundary)

simpson <- function(yval, xval) {
    stopifnot(length(xval)==length(yval))
    stopifnot(xval==sort(x=xval) || xval==rev(sort(x=xval)))
    if(length(xval)%%2==1) return(simpson_odd(yval, xval))
    if(length(xval)%%2==0) {
        n <- length(xval)
        sum1 <- simpson_odd(yval[1:(n-1)], xval[1:(n-1)]) + (xval[n]-xval[n-1])*(yval[n]+yval[n-1])/2
        sum2 <- simpson_odd(yval[2:n], xval[2:n]) + (xval[2]-xval[1])*(yval[2]+yval[1])/2
        return((sum1+sum2)/2)
    }
}

#~ interpolate_to_endpoint <- function(yval, xval, boundary) {
#~     stopifnot(length(xval)==length(yval))
#~     stopifnot(xval==sort(x=xval) || xval==rev(sort(x=xval)))
#~     if(max(xval) > boundary) {
#~         ## linearly interpolate last two points
#~     } else if (max(xval) <= boundary) {
#~         ## first determine matching interval, then 
#~     }
#~     return yb
#~ }
#~ }
