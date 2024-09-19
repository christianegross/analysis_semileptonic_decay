library("splines")

integratesplinesfromcoefficients <- function(coefficients, knots) {
    stopifnot(length(knots) <= length(coefficients[, 1]) + 1)
    sum <- 0
    for (i in seq_along(knots[1:(length(knots)-1)])) {
        xval1 <- knots[i]
        xval2 <- knots[i+1]
        part <- 0
        for (j in seq_along(coefficients[i, ])) {
            coeff <- coefficients[i, j]
            part <- part + coefficients[i, j] / j * ((xval2-xval1)^j)
        }
        sum <- sum + part
#~         print(paste(knots[i], part))
               
    }
    return (sum)
}


splineintegral <- function(yval, xval, continue=FALSE, higherlimit=0, replacelower=FALSE, lowerlimit=0, replaceindex=0) {
    stopifnot(!(continue && replacelower))
    if(any(is.na(yval))) return(NA)
    spline <- interpSpline(xval, yval)
    if(!(continue || replacelower)) {
#~         print("normal")
        return(integratesplinesfromcoefficients(spline$coefficients, spline$knots))
    }
    if(continue) {
#~         print("continue")
        stopifnot(higherlimit > max(spline$knots))
        return(integratesplinesfromcoefficients(spline$coefficients, c(spline$knots, higherlimit)))
    }
    if(replacelower) {
#~         print("replacelower")
        stopifnot(replaceindex > 0 && replaceindex <= length(xval))
        knots <- spline$knots[1:replaceindex]
        knots[replaceindex] <- lowerlimit
        return(integratesplinesfromcoefficients(spline$coefficients[1:replaceindex, ], knots))
    }
}
