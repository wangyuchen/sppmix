#' @export
selectMix <- function(pp, Ms, L = 5000, LL = 100, burnin = 1000,
                      truncate = FALSE) {
  n <- pp$n
  pattern <- cbind(pp$x,pp$y)
  lognfac <- log(sqrt(2*pi*n)) + n*log(n) - n
  mmax <- length(Ms)
  AIC <- rep(0, mmax)
  BIC <- rep(0, mmax)
  ICLC <- rep(0, mmax)
  marginal <- rep(0, mmax)
  for (m in 1:mmax) {
   post_real <- est_mix_damcmc(pp, m = Ms[m], L = L, LL = LL,
                               truncate = truncate)
   zs <- GetAvgLabelsDiscrete2Multinomial_sppmix(
     post_real$genzs[(burnin + 1):L, ], Ms[m])
   zsn0 <- zs[zs!=0]
   entropy <- sum(-zsn0*log(zsn0))
   real <- FixLS_da(post_real, pp$window)
   marginal[m] <- normmix_marginal(real)
   post_est <- get_post(real)
   mlambda <- post_est$mean_lambda
   mix <- vector("list", m)
   for(i in 1:Ms[m]){
     mix[[i]]=list(p = post_est$post_normmix$ps[i],
                   mu = post_est$post_normmix$mus[[i]],
                   sigma = post_est$post_normmix$sigmas[[i]]);
   }
    den <-densNormMix_atxy_sppmix(pattern, mix)
   loglikelihood <- -lognfac + n*log(mlambda) - mlambda +
     sum(log(den))
   r <- 1 + 6 * Ms[m]
   AIC[m] <- 2 * r - 2 * loglikelihood
   BIC[m] <- r * log(n) - 2 * loglikelihood
   ICLC[m] <- r * log(n) - 2 * loglikelihood + 2 * entropy

  }
  index <- as.matrix(expand.grid(1:mmax,1:mmax))
   bayes_factor <- apply(index,1,function(x) marginal[x[2]]/marginal[x[1]])
   bayes_factor <- matrix(bayes_factor, mmax, mmax,
                          dimnames = list(paste("denominator",Ms),
                                          paste("numerator",Ms)))
  RVAL <- list(AIC = AIC,
               BIC = BIC,
               ICLC = ICLC,
               BayesFactor = bayes_factor)
  return(RVAL)
}
