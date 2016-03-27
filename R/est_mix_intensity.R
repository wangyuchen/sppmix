#' Estimate mixture model using MCMC
#'
#' These functions fits a mixture model to spatial point pattern data using
#' DAMCMC or BDMCMC.
#'
#' @param pp Point pattern object of class \code{ppp}.
#' @param m Either number of components to fit in data augmentation MCMC or
#' maximum number of component in Birth Death MCMC.
#' @param truncate logical, indicating whether truncation is used, where the
#' component density are restricted within the domain of the point pattern.
#' @param L length of MCMC chain, default to 5000.
#' @param hyper_da hyperparameters for DAMCMC, default is (3, 1, 1).
#' @rdname est_mix
#' @examples
#' fit <- est_mix_damcmc(redwood, m = 3)
#' fit
#'
#' # see vignett "Work with 2D Normal Mixtures" for more details.
#' vignette("mcmc", package = "sppmix")
#'
#' @return A list of MCMC realizations
#' \item{allgens_List}{A list of generated nornal mixture.}
#' \item{genps}{Matrix of probabilities of components in the normal
#' mixture.}
#' \item{genmus}{Matrix of means of components in the normal
#' mixture.}
#' \item{gensigmas}{Matrix of covariance matrices of components in the normal
#'  mixture.}
#' \item{genzs}{Matrix of component indicators of points.}
#' \item{genlambdas}{Matrix of lambda (average intensity).}
#' \item{data}{The original point pattern.}
#' @export
est_mix_damcmc <- function(pp, m, truncate = FALSE,
                           L = 5000, hyper_da = c(3, 1, 1)) {

  fit <- DAMCMC2d_sppmix(points = cbind(pp$x, pp$y),
                         xlims = Window(pp)$xrange, ylims = Window(pp)$yrange,
                         m = m, truncate = truncate,
                         L = L, hyperparams = hyper_da)
  fit$data <- pp
  fit$L <- L
  fit$m <- m
  class(fit) <- "damcmc_res"
  return(fit)
}


#' @export
print.damcmc_res <- function(fit) {
  cat("Normal Mixture fit by Data Augmentation MCMC \n",
      "MCMC samples:", fit$L, "\n",
      "Number of component:", fit$m, "\n")
}


#' @param lambda1 Parameter in Poisson prior, lambda = 1 by default.
#' @param lambda2 Birth rate, lambdab2 = 10 by default.
#' @param hyper_da hyperparameters for the hierarchical prior. See 'Details'.
#' @rdname est_mix
#' @details
#' Birth-death MCMC uses the same notations as the paper from Stephens,
#' M.(2000). The definition of hyperparameters can be found in Stephens's paper
#' formula (21)-(24).
#' @references Stephens, M. "Bayesian analysis of mixture models with an unknown
#' number of components—an alternative to reversible jump methods",
#' The Annals of Statistics 2000, Vol. 28, No. 1, 40–74
#'
#' @export
est_mix_bdmcmc <- function(pp, m, truncate = FALSE,
                           lambda1 = 1, lambda2 = 10, hyper = c(5,.01,3,2,1,1),
                           L = 5000) {

  fit <- BDMCMC2d_sppmix(m, data = cbind(pp$x, pp$y),
                         xlims = Window(pp)$xrange, ylims = Window(pp)$yrange,
                         truncate = truncate,
                         lamda = lambda1, lamdab = lambda2, hyper = hyper,
                         L = L)
  fit$data <- pp
  fit$L <- L
  class(fit) <- "bdmcmc_res"
  return(fit)
}


#' @export
print.bdmcmc_res <- function(fit) {
  cat("Normal Mixture fit by Birth-death MCMC\n",
      "MCMC samples:", fit$L, "\n",
      "Number of component:", sort(unique(fit$numcomp)), "\n")
}



