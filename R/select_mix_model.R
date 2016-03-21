#' Mixture Model Selection
#'
#' Finds the best number of components by computing
#' model selection criteria, including
#' AIC (Akaike Information Criterion),
#' BIC (Bayesian Information Criterion),
#' ICLC (Integrated Classification Likelihood Criterion),
#' and the Bayes factors of models with different
#' number of components against each other.
#' The entertained model is based on the posterior means
#' of MCMC runs, for several numbers of components
#' defined in the vector Ms.
#' @param runallperms Set to 0 to use an approximation to the
#' Likelihood and Entropy within the MCMC (not affected by label switching).
#' Set to 1 to use an identifiability constraint
#' in order to remove label switching from the posterior
#' means. Set to 2 to run through all permutations and use decision
#' theory (minimize squared error loss) to obtain
#' the best permutation and then compute the posterior
#' means (BEST approach, but very slow operation). The default is 0.
#' @return A list containing the following components:
#' \item{AIC}{the values of the AIC criterion}
#' \item{BIC}{the values of the BIC criterion}
#' \item{ICLC}{the values of the ICLC criterion}
#' \item{BayesFactor}{a matrix with all the Bayes Factors}
#' \item{Marginal}{the values of the marginal density}
#' \item{LogLikelihood}{the values of the LogLikelihood}
#'
#' @references McLachlan, G., and Peel, D. (2000). Finite Mixture Models. Wiley-Interscience.
#'
#' Jasra, A., Holmes, C.C. and Stephens, D. A. (2005). Markov Chain Monte Carlo Methods and the Label Switching Problem in Bayesian Mixture. Statistical Science, 20, 50-67.
#'
#' @author Sakis Micheas, Jiaxun Chen, Yuchen Wang
#' @export
#' @examples
#' # create the true mixture
#' truemix <- normmix(ps=c(.2, .6,.2), mus=list(c(0.3, 0.3), c(0.7, 0.7), c(0.5, 0.5)),
#' sigmas = list(.005*diag(2), .001*diag(2), .001*diag(2)),
#' lambda=100,win=square(1))
#' # generate the point pattern
#' genPPP <- rsppmix(truemix)
#' plot(genPPP)
#' # compute model selection criteria
#' ModelSel=selectMix(genPPP,1:5,runallperms=0)
#' # show info
#' ModelSel
#' ModelSel=selectMix(genPPP,1:5,runallperms=1)
#' # show info
#' ModelSel
#' ModelSel=selectMix(genPPP,1:5,runallperms=2)
#' # show info
#' ModelSel
selectMix <- function(pp, Ms, L = 10000, burnin = 1000,
                      truncate = FALSE,
                      runallperms=0)
{
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
    if(runallperms==0)
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
    {
      # Start the clock!
      ptm1 <- proc.time()
      post_real<-DAMCMC2d_sppmix(cbind(pp$x, pp$y),
       Window(pp)$xrange,Window(pp)$yrange,Ms[m],L,
       truncate,c(3,1,1))
      # Stop the clock
      ptm<-proc.time() - ptm1
      cat(paste("\nComputation time in seconds:",ptm[[1]],"\n"))
      meanlamda=mean(post_real$genlamdas[(burnin+1):L])
    if(runallperms==2)
      real<<-PostGenGetBestPerm_sppmix(post_real$allgens_List)
    else
      real<<-PostGenGetBestPermIdenConstraint_sppmix(post_real)
#    marginal[m] <-  normmix_marginal(real, burnin = burnin)
#get the marginal
    cat("Calculating the entropy and the marginal...")
    if(runallperms==2)
      den <-sapply( real$permuted_gens[(burnin+1):L], densNormMix_atxy_sppmix, atxy = pattern)
    else
      den <-sapply( real$allgens_List[(burnin+1):L], densNormMix_atxy_sppmix, atxy = pattern)
    logden <- log(den)
    sumlogden <- colSums(logden)
    marginal[m] <- sum(exp(sumlogden))

#    post_est <- get_post(real)
    if(runallperms==2)
      mix_of_postmeans<<-MakeMixtureList(
        real$permuted_gens,burnin)
    else
      mix_of_postmeans<<-MakeMixtureList(
        real$allgens_List,burnin)
    den <-densNormMix_atxy_sppmix(pattern, mix_of_postmeans)

    loglikelihood[m] <- -lognfac + n*log(meanlamda) - meanlamda +
      sum(log(den))

#    cat("\npassed")
    if(runallperms==2)
    #use best permutation on the z's
      permuted_genzs=PermuteZs_sppmix(
        post_real$genzs ,real$best_perm)
    else
      permuted_genzs=post_real$genzs;
    zs <- GetAvgLabelsDiscrete2Multinomial_sppmix(
      permuted_genzs[(burnin + 1):L, ], Ms[m])
    zsn0 <- zs[zs!=0]
    entropy <- sum(-zsn0*log(zsn0))
    cat(" Done\n")
#    marginal[m] <- normmix_marginal(permgens$permuted_gens, burnin = burnin)
#    cat("\npassed1")
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
  cat(paste("AIC minimum for",Ms[which.min(AIC)],"components"))
  cat(paste("\nBIC minimum for",Ms[which.min(BIC)],"components"))
  cat(paste("\nICLC minimum for",Ms[which.min(ICLC)],"components"))
  return(RVAL)
}
