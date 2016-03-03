#' Summary of an MCMC result
#'
#' Prints a brief summary of an MCMC realization result
#' @inheritParams plot_avgsurf
#' @author Athanasios Christou Micheas, Jiaxun Chen, Yuchen Wang
#' @export
#' @examples
#' # generate data
#' mix2 <- normmix(ps=c(.4, .6), mus=list(c(0.1, 0.1), c(0.8, 0.8)),
#' sigmas=list(.02*diag(2), .01*diag(2)))
#' pp2 <- rsppmix(100,mix2,square(1))
#' # Run Data augmentation MCMC and get posterior realizations
#' post=est_mix_damcmc(pp2,L = 5000,2,truncate = F)
#' summary the posterior results
#' summary(post)
summary.damcmc_res <- function(fit, dgt = 2) {

  m = dim(fit$genmus)[1]
  for (i in 1:m)
  {
    #true value and credible sets
    poststats = GetStats_sppmix(fit$genps[,i], alpha=0.05)
    poststats <- sapply(poststats, format, digits = dgt)

    cat("\n---------------- Component ",i,"------------\n")
    #cat(paste("Probability: true =",truemix[[i]]$p))
    cat(paste("\nProbability: posterior mean =",
              poststats$Mean,"\n"))
    cat(paste(poststats$CredibleSetConfidence,
              "% Credible Set:\n[",poststats$CredibleSet[1],
              ",",poststats$CredibleSet[2],"]" ))

    poststats = GetStats_sppmix(fit$genmus[i,1,],alpha=0.05)
    poststats <- sapply(poststats, format, digits = dgt)

    #cat(paste("\nMean vector, x-coord: true =",truemix[[i]]$mu[1]))
    cat(paste("\nMean vector, x-coord: post mean =",
              poststats$Mean,"\n"))
    cat(paste(poststats$CredibleSetConfidence,
              "% Credible Set:\n[",poststats$CredibleSet[1],
              ",",poststats$CredibleSet[2],"]" ))

    poststats = GetStats_sppmix(fit$genmus[i,2,],alpha=0.05)
    poststats <- sapply(poststats, format, digits = dgt)

    #cat(paste("\nMean vector, y-coord: true =",truemix[[i]]$mu[2]))
    cat(paste("\nMean vector, y-coord: post mean =",
              poststats$Mean,"\n"))
    cat(paste(poststats$CredibleSetConfidence,
              "% Credible Set:\n[",poststats$CredibleSet[1],
              ",",poststats$CredibleSet[2],"]" ))
  }
  cat("\n----------------Component stats done------------\n")
  #cat("\nNOTE: if you have the truth, then what is called comp
}
