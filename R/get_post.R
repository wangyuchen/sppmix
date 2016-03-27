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
#' @rdname get_post
get_post <- function(obj, ...) {
  UseMethod("get_post", obj)
}


#' @examples
#'
#' fit <- est_mix_damcmc(pp = redwood, m = 3)
#' post_intsurf <- get_post(fit, burnin = 1000)
#' @rdname get_post
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


#' @param comp number of components. The posterior will be calculated only
#' based on iterations where
#' @examples
#' fit <- est_mix_bdmcmc(pp = redwood, m = 5)
#' post_intsurf <- get_post(fit, 3, burnin = 1000)
#' @rdname get_post
#' @export
get_post.bdmcmc_res <- function(fit, comp, burnin) {
  L <- length(fit$allgens)
  if (missing(burnin)) burnin <- L / 10

  numcomp <- fit$numcomp
  numcomp[1:burnin] <- 0
  ind <- numcomp == comp

  post_ps <- colMeans(fit$genps[ind, , drop = FALSE])
  post_ps <- post_ps[seq_len(comp)]

  if (any(is.nan(post_ps))) {
    stop(paste("This BDMCMC chain does not contain any", comp,
             "component realizations after burn-in."))
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
                            fit$data$window, estimated = TRUE)
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



