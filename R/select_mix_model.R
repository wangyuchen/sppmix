#' @export
selectMix <- function(pp, Ms, L = 10000, burnin = 1000,
                      truncate = FALSE) {
  if (any(Ms < 1)) {
    stop("Operation exited. Must have at least 1 component.")
  }
  n <- pp$n
  pattern <- cbind(pp$x,pp$y)
  lognfac <- log(sqrt(2*pi*n)) + n*log(n) - n
  mmax <- length(Ms)
  AIC <- rep(0, mmax)
  BIC <- rep(0, mmax)
  ICLC <- rep(0, mmax)
  marginal <- rep(0, mmax)
  loglikelihood <- rep(0, mmax)
  for (m in 1:mmax)
  {
#   post_real<-est_mix_damcmc(pp, m = Ms[m], L = L,truncate = truncate)
    cat(paste("\n================ # Components=",Ms[m],"===============\n"))
    if(1)
    {
      # Start the clock!
      ptm1 <- proc.time()
      post_real<-DAMCMC2dExtras_sppmix(cbind(pp$x, pp$y),
          Window(pp)$xrange,Window(pp)$yrange,
                  Ms[m],L,burnin,truncate,c(3,.3,1))
      # Stop the clock
      ptm<-proc.time() - ptm1
      cat(paste("\nComputation time in seconds:",ptm[[1]]))
    marginal[m]=post_real$Marginal;
    entropy=post_real$Entropy;
    loglikelihood[m]=post_real$LogLikelihood;
    }
    else
    {#too damn slow, forget about it
    post_real<-DAMCMC2d_sppmix(cbind(pp$x, pp$y),
     Window(pp)$xrange,Window(pp)$yrange,Ms[m],L,
     truncate,c(3,1,1))
    meanlamda=mean(post_real$genlamdas[(burnin+1):L])
    real <- FixLS_da(post_real)
    marginal[m] <- normmix_marginal(real, burnin = burnin)
    post_est <- get_post(real)
    mlambda <- post_est$mean_lambda
    mix <- vector("list", m)
    for(i in 1:Ms[m]){
      mix[[i]]=list(p = post_est$post_normmix$ps[i],
                    mu = post_est$post_normmix$mus[[i]],
                    sigma = post_est$post_normmix$sigmas[[i]]);
    }
    den <-densNormMix_atxy_sppmix(pattern, mix)

#    cat("\npassed")
    loglikelihood[m] <- -lognfac + n*log(meanlamda) - meanlamda +
      sum(log(den))

    zs <- GetAvgLabelsDiscrete2Multinomial_sppmix(
      post_real$genzs[(burnin + 1):L, ], Ms[m])
    zsn0 <- zs[zs!=0]
    entropy <- sum(-zsn0*log(zsn0))
#    marginal[m] <- normmix_marginal(permgens$permuted_gens, burnin = burnin)
    cat("\npassed1")
    }
    r <-6 * Ms[m]
    AIC[m] <- 2 * r - 2 * loglikelihood[m]
    BIC[m] <- r * log(n) - 2 * loglikelihood[m]
    ICLC[m] <- r * log(n) - 2 * loglikelihood[m] + 2 * entropy
    cat(paste("\nAIC=",AIC[m]))
    cat(paste("\nBIC=",BIC[m]))
    cat(paste("\nICLC=",ICLC[m]))
  }
  index <- as.matrix(expand.grid(1:mmax,1:mmax))
  bayes_factor <- apply(index,1,function(x) marginal[x[2]]/marginal[x[1]])
  bayes_factor <- matrix(bayes_factor, mmax, mmax,
                          dimnames = list(paste("denominator",Ms),
                                          paste("numerator",Ms)))
  cat("\n===============================\n")
  RVAL <- list(AIC = AIC,
               BIC = BIC,
               ICLC = ICLC,
               BayesFactor = bayes_factor,
               Marginal=marginal,
               LogLikelihood=loglikelihood)
  cat(paste("AIC minimum for",which.min(AIC),"components"))
  cat(paste("\nBIC minimum for",which.min(BIC),"components"))
  cat(paste("\nICLC minimum for",which.min(ICLC),"components"))
  return(RVAL)
}
