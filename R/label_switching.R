#' Fix label switching
#'
#' Fix label switching by apply the best permutation or identifiability constraint
#'  to the posterior realization (MCMC result)
#' @inheritParams plot_avgsurf
#' @param xlab1 The label for x-axis
#' @param ylab1 The label for y-axis
#' @param approx Logical flag indicating whether use identifiability constraint
#' to permute all realizations.
#' @param plot_result Logical flag indicating whether plot the point pattern
#' and intensity surface after permutation. The default is FALSE.
#' @author Athanasios Christou Micheas, Jiaxun Chen, Yuchen Wang
#' @export
#' @examples
#' # generate data
#' mix2 <- normmix(ps=c(.4, .6), mus=list(c(0.3, 0.3), c(0.7, 0.7)),
#' sigmas = list(.02*diag(2), .01*diag(2)))
#' pp2 <- rsppmix(100,mix2,square(1))
#' # Run Data augmentation MCMC and get posterior realizations
#' post = est_mix_damcmc(pp2,L = 5000,2,truncate = F)
#' # get posterior mean for each parameter
#' post_mean = get_post(post)
#' # plot the estimated intensity surface
#' plot(post_mean$post_normmix, post_mean$mean_lambda, square(1))
#' # Fix label switching
#' post_fixed = FixLS_da(post, plot_result = TRUE)
FixLS_da<- function(fit, burnin = length(fit$allgens_List) / 10,
                 xlab1 = "x",ylab1 = "y", approx = FALSE, plot_result = FALSE)
{
  win <- domain(fit$data)
  m <- dim(fit$genmus)[1]
  xlims1 <- c(win$xrange)
  ylims1 <- c(win$yrange)
  L <- dim(fit$genps)[1]
  if (approx == TRUE) {
    permgens <- PostGenGetBestPermIdenConstraint_sppmix(fit)
    post_ps <- colMeans(permgens$genps[-(1:burnin), ])
    mus <- apply(permgens$genmus[, , -(1:burnin)], 1:2, mean)

    mean_mat <- function(mats) Reduce("+", mats) / length(mats)
    sigmas <- apply(permgens$gensigmas[-(1:burnin), ], 2, mean_mat)

    mean_lambda <- mean(permgens$genlamdas[(burnin + 1):L])

    post_mus <- post_sigmas <- vector("list", m)
    for (i in 1:m) {
      post_mus[[i]] <- mus[i, ]
      post_sigmas[[i]] <- matrix(sigmas[, i], 2, 2)
    }
    post_normix = normmix(post_ps, post_mus, post_sigmas, lambda = mean_lambda,
                          win = win)
  } else {
    permgens <- PostGenGetBestPerm_sppmix(fit$allgens_List)
    post_ps <- colMeans(permgens$permuted_ps[-(1:burnin), ])
    mus <- apply(permgens$permuted_mus[, , -(1:burnin)], 1:2, mean)

    mean_mat <- function(mats) Reduce("+", mats) / length(mats)
    sigmas <- apply(permgens$permuted_sigmas[-(1:burnin), ], 2, mean_mat)

    mean_lambda <- mean(fit$genlamdas[(burnin + 1):L])

    post_mus <- post_sigmas <- vector("list", m)
    for (i in 1:m) {
      post_mus[[i]] <- mus[i, ]
      post_sigmas[[i]] <- matrix(sigmas[, i], 2, 2)
    }
    post_normix = normmix(post_ps, post_mus, post_sigmas, lambda = mean_lambda,
                          win = win)
  }
   if (plot_result == TRUE) {
     plot(fit$data, post_normix$mus)
     plot.intensity_surface(post_normix,
                            main = "Posterior mean intensity surface (permutated labels)")
   }
  if (approx == TRUE) {
    fit$allgens_List = permgens$allgens_List
    fit$genps = permgens$genps
    fit$genmus = permgens$genmus
    fit$gensigmas = permgens$gensigmas
    fit$genzs = permgens$genzs
    fit$genlamdas = permgens$genlamdas
  } else {
    fit$allgens_List = permgens$permuted_gens
    fit$genps = permgens$permuted_ps
    fit$genmus = permgens$permuted_mus
    fit$gensigmas = permgens$permuted_sigmas
    fit$genlamdas = fit$genlamdas
  }
   return(fit)
}

#' Test if posterior realizations of mu have label switching
#'
#' Test if there is a label switching of the posterior realizations of mu by
#' testing if the mean of each component changed dramatically during MCMC.
#'
#'@inheritParams plot_avgsurf
#' @author Athanasios (Sakis) Micheas
#' @export
test_labswitch<- function(fit, burnin = fit$L/10 ) {
  m <- fit$m
  genmus <- drop_realization(fit, burnin)$genmus
  cat("\nChecking for label switching...\n")
  for (i in 1:m) {
    if(Check4LabelSwitching_sppmix(genmus[i, 1, ])) {
      cat("Label switching present. \nPermute the labels to get a better fit,
          \nor obtain the average of the surfaces\n ")
      return(TRUE)
    }
    if(Check4LabelSwitching_sppmix(genmus[i, 2, ])) {
      cat("Label switching present. \nPermute the labels to get a better fit,
          \nor obtain the average of the surfaces\n ")
      return(TRUE)
    }
  }
  cat("No Label switching detected")
  return(invisible(FALSE))
}

