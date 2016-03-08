#' Monte Carlo goodness of fit test for normal mixture model
#'
#' Performs a test of goodness-of-fit for a given IPPP with mixture intensity
#' surface. Monte Carlo method is used to compare an observed point pattern with
#' a proposed normal mixture model.
#'
#' @param pp point pattern object of class \code{ppp}.
#' @param L length of MCMC chain, default to 5000.
#' @param burnin number of burn-in iterations.
#' @param intsurf Object of class intensity_surface.
#' @param alpha Significant level for the goodness-of-fit test.
#' @param truncate Whether to truncate the points (or mass) outside the domain,
#' default to TRUE.
#' @export
#' @examples
#' # Generate intensity surface
#' intsurf1 <- normmix(ps = c(.3, .7),
#'                     mus = list(c(0.2, 0.2), c(.8, .8)),
#'                     sigmas = list(.01*diag(2), .01*diag(2)),
#'                     lambda = 100,
#'                     win = square(1))
#' # Generate point pattern by using intsurf1
#' pp1 <- rsppmix(intsurf1)
#' # Apply goodness-of-fit test
#' mc_gof(pp1, intsurf1, 0.05)
#' # Generate a different intensity surface
#' intsurf2 <- normmix(ps = c(.5, .5),
#'                     mus = list(c(0.2, 0.8), c(.8, .2)),
#'                     sigmas = list(.01*diag(2), .01*diag(2)),
#'                     lambda = 100,
#'                     win = square(1))
#' # Apply goodness-of-fit test
#' mc_gof(pp1, intsurf2, 0.05)

mc_gof <- function(pp, intsurf, alpha, L = 5000, burnin = L/10,
                   truncate = TRUE) {
  # check first
  if (L < 1000) {
    stop("L need to be larger than 1000.")
  }
  n <- pp$n
  if (n == 0) {
    stop("No points.")
  }
  m <- intsurf$m
  if (m == 0) {
    stop("No components in the mixture.")
  }
  win <- intsurf$window
  # get posterior realizations
  post <- est_mix_damcmc(pp = pp, m = m, truncate = truncate, L = L)
  # post_mix <- get_post(post)
  T_max <- rep(0, (L - burnin))
  T_mean <- rep(0, (L - burnin))
  T_var <- rep(0, (L - burnin))
  # get posterior predictive sample
  for (i in 1:(L-burnin)) {
    lambda <- post$genlamdas[i+burnin]
    ps <- post$genps[(i+burnin), ]
    mus <- split(post$genmus[, , (i+burnin)], 1:m)
    sigma <- post$gensigmas[(i+burnin), ]
    mix_real <- normmix(ps, mus, sigma, lambda = lambda, win = win)
    ow <- options("warn")
    options(warn = -1)
    pp_pred <- rsppmix(mix_real, truncate = truncate)
    options(ow)
    # compute the test statistics for the predicted pattern
    summean <- 0
    summax <- 0
    sumvar <- 0
    ind <- unique(pp_pred$comp)
    if (length(ind) != m) {
      mispp <- setdiff(1:m, ind)
      message(paste("No point generated from component ", mispp))
    }
    for (k in ind) {
      ppk <- pp_pred[pp_pred$comp == k]
      distk <- sqrt(rowSums((cbind(ppk$x, ppk$y) - mus[[k]])^2))
      summean <- summean + mean(distk)
      summax <- summax + max(distk)
      sumvar <- sumvar + var(distk)
    } # end loop of k
    T_mean[i] <- summean/m
    T_max[i] <- summax/m
    T_var[i] <- sumvar/m
  } # end look of i
  sortT_mean <- sort(T_mean)
  sortT_max <- sort(T_max)
  sortT_var <- sort(T_var)
  # compute test statistcs for the given pattern
  summean <- 0
  summax <- 0
  sumvar <- 0
  probz <- GetAvgLabelsDiscrete2Multinomial_sppmix(post$genzs[(burnin + 1):L, ],
                                                   m)
    for (j in 1:m) {
    distj <- probz[, j]*sqrt(rowSums((cbind(pp$x, pp$y) - intsurf$mus[[j]])^2))
    summean <- summean + mean(distj)
    summax <- summax + max(distj)
    sumvar <- sumvar + var(distj)
  } # end loop of j
  T_meanTS <- summean/m
  T_maxTS <- summax/m
  T_varTS <- sumvar/m
  # test for T_mean
  T_meanalpha <- sortT_mean[floor((1 - alpha)*(L - burnin))]
  # Ts_mean <- T_meanTS > T_meanalpha
  # if (Ts_mean == TRUE) {
  #   result_mean <- "One-sided GOF using T_meanTS, Reject Null.\\n
  #   Mixture Model DOES NOT FIT WELL"
  # } else {
  #   result_mean <- "One-sided GOF using T_meanTS, Cannot Reject Null.\\n
  #   Mixture Model FITS WELL"
  # }
  p_mean <- mean(T_meanTS < T_mean)

  # test for T_max
  T_maxalpha <- sortT_max[floor((1 - alpha)*(L - burnin))]
  # Ts_max <- T_maxTS > T_maxalpha
  # if (Ts_max == TRUE) {
  #   result_max <- "One-sided GOF using T_maxTS, Reject Null.\\n
  #   Mixture Model DOES NOT FIT WELL"
  # } else {
  #   result_max <- "One-sided GOF using T_maxTS, Cannot Reject Null.\\n
  #   Mixture Model FITS WELL"
  # }
  p_max <- mean(T_maxTS < T_max)

  # test for T_var
  T_varalpha <- sortT_var[floor((1 - alpha)*(L - burnin))]
  # Ts_var <- T_varTS > T_varalpha
  # if (Ts_var == TRUE) {
  #   result_var <- "One-sided GOF using T_varTS, Reject Null.\\n
  #   Mixture Model DOES NOT FIT WELL"
  # } else {
  #   result_var <- "One-sided GOF using T_varTS, Cannot Reject Null.\\n
  #   Mixture Model FITS WELL"
  # }
  p_var <- mean(T_varTS < T_var)
  RVAL <- list()
  RVAL[[1]] <- list(statistic=c(T_meanTS = T_meanTS, T_alpha = T_meanalpha),
               p.value = p_mean,method = "One-sided GOF using T_meanTS",
               alternative="Mixture Model DOES NOT FIT WELL")
  class(RVAL[[1]]) <- "htest"
  RVAL[[2]] <- list(statistic=c(T_maxTS = T_maxTS, T_alpha = T_maxalpha),
                    p.value = p_max,method = "One-sided GOF using T_maxTS",
                    alternative="Mixture Model DOES NOT FIT WELL")
  class(RVAL[[2]]) <- "htest"
  RVAL[[3]] <- list(statistic=c(T_varTS = T_varTS, T_alpha = T_varalpha),
                    p.value = p_var,method = "One-sided GOF using T_varTS",
                    alternative="Mixture Model DOES NOT FIT WELL")
  class(RVAL[[3]]) <- "htest"
  return(RVAL)
}
