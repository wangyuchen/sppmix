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
#' @param LL additional fitting parameter
#'
#' @rdname est_mix
#' @examples
#' # see vignett "Work with 2D Normal Mixtures" for details.
#'
#' @export
est_mix_damcmc <- function(pp, m, truncate = FALSE,
                           L = 5000, LL = 100) {

  fit <- DAMCMC2d_sppmix(points = cbind(pp$x, pp$y),
                         xlims = Window(pp)$xrange, ylims = Window(pp)$yrange,
                         m = m, truncate = truncate,
                         L = L, LL = LL)
  fit$data <- pp
  class(fit) <- "damcmc_res"
  return(invisible(fit))
}

#' Calculate posterior mixture from estimated mixture
#'
#' @param fit object from \code{est_mix_damcmc} or \code{est_mix_bdmcmc}.
#' @param burnin number of burn-in iterations.
#'
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

  return(list(post_normmix = normmix(post_ps, post_mus, post_sigmas),
              mean_lambda = mean_lambda))
}

#' @param lambda parameter for BDMCMC
#' @param lambdab parameter for BDMCMC
#' @param hyper parameter for BDMCMC
#' @rdname est_mix
#' @export
est_mix_bdmcmc <- function(pp, m, truncate = FALSE,
                           lambda, lambdab, hyper,
                           L = 5000, LL = 100) {

  fit <- BDMCMC2d_sppmix(m, data = cbind(pp$x, pp$y),
                         xlims = Window(pp)$xrange, ylims = Window(pp)$yrange,
                         truncate = truncate,
                         lamda = lambda, lamdab = lambdab, hyper = hyper,
                         L = L, LL = LL)
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
  ind <- fit$numcomp == comp

  comp_rlz <- GetBDCompRealiz_sppmix(fit$allgens_List[-(1:burnin)],
                                     fit$genlamdas[-(1:burnin)],
                                     fit$numcomp[-(1:burnin)], comp)


  post_ps <- colMeans(fit$genps[ind, , drop = FALSE])
  post_ps <- post_ps[seq_len(comp)]

  mus <- apply(fit$genmus[, , ind, drop = FALSE], 1:2, mean)
  mus <- mus[seq_len(comp), ]

  mean_mat <- function(mats) Reduce("+", mats) / length(mats)
  sigmas <- apply(fit$gensigmas[ind, , drop = FALSE], 2, mean_mat)
  sigmas <- sigmas[seq_len(comp)]

  mean_lambda <- mean(fit$genlamdas[ind])

  post_mus <- vector("list", comp)
  for (i in 1:comp) {
    post_mus[[i]] <- mus[i, ]
  }

  return(list(post_normmix = normmix(post_ps, post_mus, sigmas),
              mean_lambda = mean_lambda))
}










