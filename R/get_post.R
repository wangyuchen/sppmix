#' Get Posterior normal mixture
#'
#' Calculate posterior mixture after DAMCMC or BDMCMC fit. The posterior
#' mixture is calculated by using posterior mean of parameters. For birth-death
#' MCMC, the number of component need to be specified.
#'
#'
#' @param fit object from \code{est_mix_damcmc} or \code{est_mix_bdmcmc}.
#' @param burnin number of burn-in iterations. By default, it's 1/10 of chain
#' length.
#'
#' @return An object of class \code{intensity_surface}.
#'
#' @export
#' @rdname get_post
get_post <- function(obj, ...) {
  UseMethod("get_post", obj, ...)
}


#' @examples
#'
#' fit <- est_mix_damcmc(pp = redwood, m = 3)
#' post_intsurf <- get_post(fit, burnin = 1000)
#' @rdname get_post
#' @export
get_post.damcmc_res <- function(fit, burnin = fit$L / 10) {
  fit_burnined <- drop_realization(fit, burnin)

  post_ps <- colMeans(fit_burnined$genps[, 1:fit$m])
  mus <- apply(fit_burnined$genmus[1:fit$m, , ], 1:2, mean)

  sigmas <- apply(fit_burnined$gensigmas[, 1:fit$m], 2,
                  function(mats) Reduce(`+`, mats) / length(mats))

  mean_lambda <- mean(fit_burnined$genlamdas)

  post_mus <- post_sigmas <- vector("list", fit$m)
  for (i in seq_along(post_mus)) {
    post_mus[[i]] <- mus[i, ]
    post_sigmas[[i]] <- matrix(sigmas[, i], 2, 2)
  }
  normmix(post_ps, post_mus, post_sigmas, mean_lambda,
          fit$data$window, estimated = TRUE)
}


#' @param comp number of components. The posterior will be calculated only
#' based on iterations where
#' @examples
#' fit <- est_mix_bdmcmc(pp = redwood, m = 5)
#' post_intsurf <- get_post(fit, 3, burnin = 1000)
#' @rdname get_post
#' @export
get_post.bdmcmc_res <- function(fit, comp, burnin = fit$L / 10) {
  fit_burnined <- drop_realization(fit, burnin)

  if (!(comp %in% unique(fit_burnined$numcomp))) {
    stop(paste("This BDMCMC chain does not contain any", comp,
               "component realizations after burn-in."))
  }

  new_chain <- drop_realization(fit_burnined, fit_burnined$numcomp != comp)
  new_chain$m <- comp

  get_post.damcmc_res(new_chain, burnin = 0)

}


drop_realization <- function(fit, drop) {
  # burn-in: subsetting vector (integer or logical) for the realizations
  # to be dropped
  if (is.numeric(drop)) {
    if (length(drop) == 1) keep <- seq_len(fit$L) > drop else keep <- -drop
  }

  if (is.logical(drop)) keep <- !drop

  fit$allgens_List <- fit$allgens_List[keep]
  fit$genps <- fit$genps[keep, , drop = FALSE]
  fit$genmus <- fit$genmus[, , keep, drop = FALSE]
  fit$gensigmas <- fit$gensigmas[keep, , drop = FALSE]
  fit$genzs <- fit$genzs[keep, , drop = FALSE]
  fit$genlamdas <- fit$genlamdas[keep, , drop = FALSE]
  fit$ApproxCompMass <- fit$ApproxCompMass[keep, , drop = FALSE]

  fit$L <- length(fit$allgens_List)
  if (class(fit) == "bdmcmc_res")
    fit$numcomp <- fit$numcomp[keep, , drop = FALSE]

  fit
}




