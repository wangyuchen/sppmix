% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kstest2d.R
\name{kstest2d}
\alias{kstest2d}
\title{Nonparametric Goodness-of-fit test for two-dimesional point patterns.}
\usage{
kstest2d(x1, x2)
}
\arguments{
\item{x1, x2}{Objects of class \code{\link[spatstat]{ppp}}.}
}
\value{
A list with class "htest" containing the following components:
\item{statistic}{Value of the KS statistic}
\item{p.value}{The p-value of the test}
\item{alternative}{A character string describing the alternative hypothesis}
}
\description{
This function performs a two-dimensional Kolmogorov-Smirnov goodness-of-fit
test of two point patterns.
}
\examples{
# genrating two point patterns
mix1 <- rnormmix(3, sig0 = .01, df = 5, xlim = c(0, 5), ylim = c(0, 5))
mix2 <- rnormmix(8, sig0 = .01, df = 10, xlim = c(0, 5), ylim = c(0, 5))
pp1 <- rsppmix(mix1, lambda = 20, win = square(5))
pp2 <- rsppmix(mix2, lambda = 20, win = square(5))

# Test for goodness of fit
kstest2d(pp1, pp2)
}
\references{
J. A. Peacock, "Two-dimensional goodness-of-fit testing in
astronomy", Monthly Notices Royal Astronomy Society 202 (1983)615-627.

Adapted from Matlab code by Dylan Muir.
}

