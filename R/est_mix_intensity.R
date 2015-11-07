#' @export
est_mix_damcmc <- function(pp, m, truncate = FALSE,
                           L = 5000, LL = 100) {

  fit <- DAMCMC2d_sppmix(data = cbind(pp$x, pp$y),
                         xlims = Window(pp)$xrange, ylims = Window(pp)$yrange,
                         m = m, truncate = truncate,
                         L = L, LL = LL)
  fit$data <- pp
  class(fit) <- "damcmc_res"
  return(invisible(fit))
}


#' @export
get_post <- function(obj, ...) {
  UseMethod("get_post", obj)
}

#' @export
get_post.damcmc_res <- function(fit, burnin) {
  L <- length(fit$allgens)
  m <- ncol(fit$genps)
  if (missing(burnin)) burnin <- L / 10

  post_ps <- colMeans(fit$genps[-(1:burnin), ])
  mus <- apply(fit$genmus[, , -(1:burnin)], 1:2, mean)

  mean_mat <- function(mats) Reduce("+", mats) / length(mats)
  sigmas <- apply(fit$gensigmas[-(1:burnin), ], 2, mean_mat)

  mean_lambda <- mean(fit$genlamdas)

  post_mus <- post_sigmas <- vector("list", m)
  for (i in 1:m) {
    post_mus[[i]] <- mus[i, ]
    post_sigmas[[i]] <- matrix(sigmas[, i], 2, 2)
  }

  return(list(post_normmix = normmix(post_ps, post_mus, post_sigmas),
              mean_lambda = mean_lambda))
}


#' @export
est_mix_bdmcmc <- function(pp, max_comp, truncate = FALSE,
                           lambda, lambdab, hyper,
                           L = 5000, LL = 100) {

  fit <- BDMCMC2d_sppmix(max_comp, data = cbind(pp$x, pp$y),
                         xlims = Window(pp)$xrange, ylims = Window(pp)$yrange,
                         truncate = truncate,
                         lamda = lambda, lamdab = lambdab, hyper = hyper,
                         L = L, LL = LL)
  class(fit) <- "bdmcmc_res"
  return(invisible(fit))
}
