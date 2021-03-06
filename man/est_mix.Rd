% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est_mix.R
\name{est_mix}
\alias{est_mix}
\title{Estimate mixture models using MCMC.}
\usage{
est_mix(pattern, m, method = c("DAMCMC"), truncate = TRUE, chains = 4,
  iter = 5000, burnin = iter/5, thin = 1, ...)
}
\arguments{
\item{pattern}{Point pattern object of class \code{ppp}.}

\item{m}{Either number of components to fit in data augmentation MCMC or
maximum number of component in Birth Death MCMC.}

\item{method}{Character specifying which method to use for fitting mixture.
Currently only "DAMCMC" is supported.}

\item{truncate}{logical, indicating whether truncation is used, where the
component density are restricted within the domain of the point pattern.}

\item{chains}{A positive integer specifying number of chains; defaults to 4.}

\item{iter}{Number of iterations for each chain (including burnin), default
to 5000.}

\item{burnin}{Number of burn-in iterations, default to 1/5 of chain length.}

\item{thin}{Positive integer for thinning parameter, default to 1.}

\item{...}{Further parameters to be passed to specific methods.}
}
\description{
This functions fits a mixture model to spatial point pattern data using
DAMCMC or BDMCMC.
}
\examples{
fit <- est_mix_damcmc(redwood, m = 3)
fit

# see vignett "Work with 2D Normal Mixtures" for more details.
vignette("mcmc", package = "sppmix")

}

