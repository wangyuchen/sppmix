#' Get Posterior normal mixture.
#'
#' Calculate posterior normal mixture of a DAMCMC or BDMCMC fit. The posterior
#' mixture is calculated by using posterior mean of parameters. For birth-death
#' MCMC, the number of component should be specified, and all realization with
#' that number of components are gathered to calculate posterior mixture.
#'
#' @param fit Object of class \code{damcmc_res} or \code{bdmcmc_res}.
#' @param burnin Number of burn-in iterations. By default, it's 1/10 of chain
#' length.
#'
#' @return An object of class \code{intensity_surface}.
#'
#' @export
#' @rdname get_post
get_post <- function(obj, ...) {
  UseMethod("get_post", obj)
}


#' @examples
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


#' @param num_comp Number of components. The posterior will be calculated only
#' based on iterations have that many number of components.
#' @examples
#'
#' fit <- est_mix_bdmcmc(pp = redwood, m = 5)
#' post_intsurf <- get_post(fit, num_comp = 4, burnin = 1000)
#' @rdname get_post
#' @export
get_post.bdmcmc_res <- function(fit, num_comp, burnin = fit$L / 10) {
  fit_burnined <- drop_realization(fit, burnin)

  if (!(num_comp %in% unique(fit_burnined$numcomp))) {
    stop(paste0("This BDMCMC chain does not contain any ", num_comp,
               "-component realizations after burn-in."))
  }

  new_chain <- drop_realization(fit_burnined, fit_burnined$numcomp != num_comp)
  new_chain$m <- num_comp

  get_post.damcmc_res(new_chain, burnin = 0)
}

#' Drop a subset of MCMC realization
#'
#' Drop realizations from a MCMC chain.
#'
#' @inheritParams get_post
#' @param drop If one integer is provided, it will drop the first 1:drop
#' realizations. If an integer vector is provided, it will drop these
#' iterations. If a logical vector is provided (with the same length as chain
#' length of \code{fit}), it will be used for subsetting directly.
#'
#' @examples
#' fit <- est_mix_bdmcmc(redwood, m = 5)
#' fit
#' drop_realization(fit, 500)
#' drop_realization(fit, fit$numcompo != 5)
#'
drop_realization <- function(fit, drop) {
  # Generalized function for burn-in or drop realizations
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
  if (class(fit) == "bdmcmc_res") {
    fit$numcomp <- fit$numcomp[keep, , drop = FALSE]
  }

  fit
}




