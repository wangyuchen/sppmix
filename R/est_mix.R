#' Estimate mixture models using MCMC.
#'
#' This functions fits a mixture model to spatial point pattern data using
#' DAMCMC or BDMCMC.
#'
#' @param pattern Point pattern object of class \code{ppp}.
#' @param m Either number of components to fit in data augmentation MCMC or
#' maximum number of component in Birth Death MCMC.
#' @param method Character specifying which method to use for fitting mixture.
#' Currently only "DAMCMC" is supported.
#' @param truncate logical, indicating whether truncation is used, where the
#' component density are restricted within the domain of the point pattern.
#' @param chains A positive integer specifying number of chains; defaults to 4.
#' @param iter Number of iterations for each chain (including burnin), default
#' to 5000.
#' @param burnin Number of burn-in iterations, default to 1/5 of chain length.
#' @param thin Positive integer for thinning parameter, default to 1.
#' @param ... Further parameters to be passed to specific methods.
#' @examples
#' fit <- est_mix_damcmc(redwood, m = 3)
#' fit
#'
#' # see vignett "Work with 2D Normal Mixtures" for more details.
#' vignette("mcmc", package = "sppmix")
#'
#' @export
est_mix <- function(pattern, m, method = c("DAMCMC"), truncate = TRUE,
                    chains = 4, iter = 5000, burnin = iter / 5,
                    thin = 1, ...) {

}
