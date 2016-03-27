#' Estimate mixture model using MCMC
#'
#' These functions fits a mixture model to spatial point pattern data using
#' DAMCMC or BDMCMC.
#'
#' @param pp point pattern object of class \code{ppp}.
#' @param m either number of components to fit in data augmentation MCMC or
#' maximum number of component in Birth Death MCMC.
#' @param truncate logical, indicating whether truncation is used, where the
#' component density are restricted within the domain of the point pattern.
#' @param L length of MCMC chain, default to 5000.
#' @param hyper_da hyperparameters for DAMCMC, default is (3, 1, 1).
#' @rdname est_mix
#' @examples
#' # see vignett "Work with 2D Normal Mixtures" for details.
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
#' \item{genlambdas}{Matrix of lambda(average intensity).}
#' \item{data}{The original point pattern.}
#' @export
est_mix_damcmc <- function(pp, m, truncate = FALSE,
                           L = 5000, hyper_da = c(3, 1, 1)) {

  fit <- DAMCMC2d_sppmix(points = cbind(pp$x, pp$y),
                         xlims = Window(pp)$xrange, ylims = Window(pp)$yrange,
                         m = m, truncate = truncate,
                         L = L, hyperparams = hyper_da)
  fit$data <- pp
  class(fit) <- "damcmc_res"
  return(invisible(fit))
}

#' Calculate posterior mixture by using posterior mean of parameters
#'
#' @param fit object from \code{est_mix_damcmc} or \code{est_mix_bdmcmc}.
#' @param burnin number of burn-in iterations.
#'
#' @return An object of class \code{intensity_surface}.
#'
#' @details
#' For Birth-Death MCMC, the number of component need to be specified.
#' @export
#' @examples
#'
#' fit <- sppmix::est_mix_damcmc(pp = redwood, m = 3, truncate = FALSE,
#'                               L = 50000, LL = 100)
#' get_post(fit, burnin = 5000)
#'
get_post <- function(obj, ...) {
  UseMethod("get_post", obj)
}

#' @export
get_post.damcmc_res <- function(fit, burnin) {
  L <- length(fit$allgens)
  m <- ncol(fit$genps)
  if (missing(burnin)) burnin <- L / 10

  post_ps <- colMeans(fit$genps[-(1:burnin), , drop = FALSE])
  mus <- apply(fit$genmus[, , -(1:burnin), drop = FALSE], 1:2, mean)

  mean_mat <- function(mats) Reduce("+", mats) / length(mats)
  sigmas <- apply(fit$gensigmas[-(1:burnin), , drop = FALSE], 2, mean_mat)

  mean_lambda <- mean(fit$genlamdas)

  post_mus <- post_sigmas <- vector("list", m)
  for (i in 1:m) {
    post_mus[[i]] <- mus[i, ]
    post_sigmas[[i]] <- matrix(sigmas[, i], 2, 2)
  }
  post_intensity <- normmix(post_ps, post_mus, post_sigmas, mean_lambda,
                            fit$data$window, estimated = TRUE)
  return(post_intensity)
}

#' @param lambda1 Parameter in Poisson prior, lambda = 1 by default.
#' @param lambda2 Birth rate, lambdab2 = 10 by default.
#' @param hyper_da hyperparameters for the hierarchical prior. See 'Details'.
#' @rdname est_mix
#' @details
#' Birth-Death MCMC uses the same notations as the paper from Stephens, M.(2000).
#' The definition of hyperparameters can be found in Stephens's paper formula (21)-(24).
#' @references Stephens, M. "Bayesian analysis of mixture models with an unknown
#' number of components—an alternative to reversible jump methods", The Annals of Statistics 2000,
#' Vol. 28, No. 1, 40–74
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
  class(fit) <- "bdmcmc_res"
  return(invisible(fit))
}




#' @export
get_post.bdmcmc_res <- function(fit, comp, burnin) {
  L <- length(fit$allgens)
  if (missing(burnin)) burnin <- L / 10

  numcomp <- fit$numcomp
  numcomp[1:burnin] <- 0
  ind <- numcomp == comp

  # comp_rlz <- GetBDCompRealiz_sppmix(fit$allgens_List,
  #                                    fit$genlamdas,
  #                                    fit$numcomp, comp)


  post_ps <- colMeans(fit$genps[ind, , drop = FALSE])
  post_ps <- post_ps[seq_len(comp)]

  if (any(is.nan(post_ps))) {
    stop(paste("This BDMCMC does not contain any ", comp, " components
               realizations"))
  }

  mus <- apply(fit$genmus[, , ind, drop = FALSE], 1:2, mean)
  mus <- matrix(mus[seq_len(comp), ], comp, 2)

  mean_mat <- function(mats) Reduce("+", mats) / length(mats)
  sigmas <- apply(fit$gensigmas[ind, , drop = FALSE], 2, mean_mat)
  sigmas <- sigmas[seq_len(comp)]

  mean_lambda <- mean(fit$genlamdas[ind])

  post_mus <- vector("list", comp)
  for (i in 1:comp) {
    post_mus[[i]] <- mus[i, ]
  }
  post_intensity <- normmix(post_ps, post_mus, sigmas, mean_lambda,
                                      fit$data$window)
  return(post_intensity)
}


#' @export
get_bdmcmc_rlz <- function(fit, comp) {

  L <- length(fit$allgens)

  res <- list()
  ind <- fit$numcomp == comp

  res$allgens_List <- fit$allgens_List[ind]
  res$genps <- fit$genps[ind, seq_len(comp), drop = FALSE]
  res$genmus <- fit$genmus[seq_len(comp), , ind, drop = FALSE]

  res$gensigmas <- fit$gensigmas[ind, seq_len(comp), drop = FALSE]

  res$genlamdas <- fit$genlamdas[ind, 1, drop = FALSE]
  res$genzs <- fit$genzs[ind, , drop = FALSE]
  res$data <- fit$data

  class(res) <- "damcmc_res"
  return(res)
}








