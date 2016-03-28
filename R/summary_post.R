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
summary.damcmc_res <- function(fit, burnin = fit$L / 10, dgt = 2) {
  fit <- drop_realization(fit, burnin)
  m = fit$m
  for (i in 1:m) {
    #true value and credible sets
    poststats = GetStats_sppmix(fit$genps[,i], alpha=0.05)
    poststats <- sapply(poststats, format, digits = dgt)

    cat("---------------- Component ", i , "----------------\n")
    cat("Probability: posterior mean =",
              poststats$Mean,"\n")
    cat(poststats$CredibleSetConfidence,
              "% Credible Set: [",poststats$CredibleSet[1],
              ",",poststats$CredibleSet[2],"]\n" )

    poststats = GetStats_sppmix(fit$genmus[i,1,],alpha=0.05)
    poststats <- sapply(poststats, format, digits = dgt)

    cat("Mean vector, x-coord: post mean =",
              poststats$Mean,"\n")
    cat(poststats$CredibleSetConfidence,
              "% Credible Set: [",poststats$CredibleSet[1],
              ",",poststats$CredibleSet[2],"]\n" )

    poststats = GetStats_sppmix(fit$genmus[i,2,],alpha=0.05)
    poststats <- sapply(poststats, format, digits = dgt)

    cat("Mean vector, y-coord: post mean =",
              poststats$Mean,"\n")
    cat(poststats$CredibleSetConfidence,
              "% Credible Set: [",poststats$CredibleSet[1],
              ",",poststats$CredibleSet[2],"]\n" )
  }
}
