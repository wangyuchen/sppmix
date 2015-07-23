#' Two-sample Two-dimesional Kolmogorov-Smirnov test
#'
#' Performs a test of goodness-of-fit of two point patterns.  The emperical
#' distribution of two point patterns are compared using the Kolmogorov-Smirnov
#' test.
#'
#' @param x1,x2 Objects of class \code{\link[spatstat]{ppp}}.
#'
#' @return A list with class "htest" containing the following components:
#' \item{statistic}{the value of the KS statistic}
#' \item{p.value}{the p-value of the test}
#' \item{alternative}{a character string describing the alternative hypothesis}
#'
#' @references J. A. Peacock, "Two-dimensional goodness-of-fit testingin
#' astronomy", Monthly Notices Royal Astronomy Society 202 (1983)615-627.
#'
#' Adapted from Matlab code by Dylan Muir.
#' @export
#' @examples
#' # genrating two point patterns
#' if(require(spatstat)){
#'   mix1 <- rnormmix(3, .01, 5, square(5))
#'   mix2 <- rnormmix(8, .01, 10, square(5))
#'   pp1=rsppmix(20,mix1,square(5))
#'   pp2=rsppmix(20,mix2,square(5))
#' }
#' # Test for goodness of fit
#' kstest2d(pp1,pp2)
kstest2d <- function(x1, x2) {
  ecdf2d <- function(x,edge){
    count=cbind(as.numeric(x$x>=edge[1] & x$y>=edge[2]),
                as.numeric(x$x<=edge[1] & x$y>=edge[2]),
                as.numeric(x$x<=edge[1] & x$y<=edge[2]),
                as.numeric(x$x>=edge[1] & x$y<=edge[2]))
    return(count)
  }

  if((!spatstat::is.ppp(x1)) | (!spatstat::is.ppp(x2))){
    stop("Point pattern must be of class ppp.")
  }
  n1 <- length(x1$x)
  n2 <- length(x2$x)
  pp1 <- cbind(x1$x,x1$y)
  pp2 <- cbind(x2$x,x2$y)
  KSstat=-Inf
  for (i in 1:(n1+n2)){
    if(i <= n1) {
      edge=pp1[i,]
    } else {
      edge=pp2[i-n1,]
    }
    vfcdf1=apply(ecdf2d(x1,edge), 2, sum)/n1
    vfcdf2=apply(ecdf2d(x2,edge), 2, sum)/n2
    vfdiff=abs(vfcdf2-vfcdf1)
    fksts=max(vfdiff)
    if (fksts>KSstat) KSstat=fksts
  }

  n <- n1*n2/(n1+n2)
  Zn <- sqrt(n) * KSstat
  Zinf <- Zn/(1-0.53 * n^(-0.9))
  pValue=2*exp(-2 * (Zinf - 0.5)^2)

  RVAL <- list(statistic=c(KSstat = KSstat),p.value = pValue,method =
                 "Two sample 2D Kolmogorov-Smirnov Test",data.name =
                 paste(sep=", ",deparse(substitute(x1)),deparse(substitute(x2))
                       ),alternative="No goodness of fit")
  class(RVAL) <- "htest"
  return(RVAL)
}

