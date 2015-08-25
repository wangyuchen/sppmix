#' @export
selectMix <- function(pattern, win, Ms, L = 1000,
                        burnin = 200, truncate = TRUE){
  if (truncate==TRUE) {
    pattern <- pattern[spatstat::inside.owin(pattern$x, pattern$y, win)]
  }

  pp <- as.data.frame(pattern)
  n <- npoints(pattern)

  if (L <= burnin) {
    stop("wrong L or burnin")
  }
  window_area <- diff(win$x) * diff(win$y)
  lambdahat <- n / window_area
  mmax <- length(Ms)
  AIC <- rep(0, mmax)
  BIC <- rep(0, mmax)
  marginal <- rep(0, mmax)
  models <- list()
  for (m in 1:mmax) {
   mix <- est_mix_intensity(pattern, win, m = Ms[m], marginal = TRUE)
   models <- append(models, mix)
   loglikelihood <- n*log(mix$lambda) - mix$lambda
   den <- dnormmix(mix$post_mix,win = win)$v
   loglikelihood <- loglikelihood + sum(log(den))
   r <- 1 + 6 * Ms[m]
   AIC[m] <- 2 * r - 2 * loglikelihood
   BIC[m] <- r * log(n) - 2 * loglikelihood
   marginal[m] <- mix$marginal
  }
  index <- as.matrix(expand.grid(1:mmax,1:mmax))
  bayes_factor <- apply(index,1,function(x) marginal[x[2]]/marginal[x[1]])
  bayes_factor <- matrix(bayes_factor, mmax, mmax,
                         dimnames = list(paste("denominator",Ms),
                                         paste("numerator",Ms)))
  RVAL <- list(AIC = AIC,
               BIC = BIC,
               bayes_factor = bayes_factor,
               models = models)
}
