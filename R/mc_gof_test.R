#' Monte Carlo goodness of fit test for normal mixture model
#'
#' Performs a test of goodness-of-fit for a given IPPP with mixture intensity
#' surface. Monte Carlo method is used to compare an observed point pattern with
#' a proposed normal mixture model.
#'
#' @param pp point pattern object of class \code{ppp}.
#' @param L length of MCMC chain, default to 10000.
#' @param burnin number of burn-in iterations.
#' @param intsurf Object of class intensity_surface.
#' @param alpha Significant level for the goodness-of-fit test.
#' @param truncate Whether to truncate the points (or mass) outside the domain,
#' default to TRUE.
#' @details The test statistic is the average distance between the given point patten
#' and each component mean. MCMC realizations are used to simulate the distribution
#' of the test statistic. To conduct the test, need to compare the test statistic
#' with the percentile of the simulated distribution.
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

mc_gof <- function(pp, intsurf, alpha, L = 10000, burnin = L/10,
                   truncate = FALSE) {
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
  zeronj <- 0
  win <- intsurf$window
  # get posterior realizations
  post <- est_mix_damcmc(pp = pp, m = m, truncate = truncate, L = L)
  # post_mix <- get_post(post)
  T_mean <- rep(0, (L - burnin))
  # get posterior predictive sample
  for (i in 1:(L-burnin)) {
    lambda <- post$genlamdas[i+burnin]
    ps <- post$genps[(i+burnin), ]/sum(post$genps[(i+burnin), ])
    mus <- split(post$genmus[, , (i+burnin)], 1:m)
    sigma <- post$gensigmas[(i+burnin), ]
    mix_real <- normmix(ps, mus, sigma, lambda = lambda, win = win)
    ow <- options("warn")
    options(warn = -1)
    pp_pred <- rsppmix(mix_real, truncate = truncate)
    options(ow)
    # compute the test statistics for the predicted pattern
    summean <- 0
    ind <- unique(pp_pred$comp)
    if (length(ind) != m) {
      zeronj <- zeronj + 1
    }
      for (k in ind) {
        ppk <- pp_pred[pp_pred$comp == k]
        distk <- sqrt(rowSums((sweep(cbind(ppk$x, ppk$y), 2, mus[[k]]))^2))
        summean <- summean + mean(distk)
      } # end loop of k
    T_mean[i] <- summean/m
  } # end look of i
  if (zeronj != 0) {
    message(paste("There are", zeronj, "cases with components have 0 points in them.\n",
                  "MC test could be reporting wrong.") )
  }
  sortT_mean <- sort(T_mean)
  # compute test statistcs for the given pattern
  summean <- 0
  probz <- GetAvgLabelsDiscrete2Multinomial_sppmix(
    post$genzs[(burnin + 1):L, ], m)
    for (j in 1:m) {
    distj <- probz[, j]*sqrt(rowSums((sweep(cbind(pp$x, pp$y), 2, intsurf$mus[[j]]))^2))
    summean <- summean + mean(distj)
  } # end loop of j
  T_meanTS <- summean/m
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
  if (p_mean < alpha) {
    result_mean <- "\nSmall p-value, Mixture Model does not fit well.\n "
  } else {
    result_mean <- " Large p-value, Mixture Model fits the data well. \n"
  }
  result <- paste("\n         Monte Carlo Goodness of Fit Test \nTest Statistic:",
                  round(T_meanTS, 4), " Critical Value:", round(T_meanalpha, 4),
                  "\nNull Hypothesis: Mixture Model fits the pattern well.",
                  "\nAlternative: Mixture Model is not appropriate. \np-value:",
                  round(p_mean, 4), result_mean)
  cat(result)
}
