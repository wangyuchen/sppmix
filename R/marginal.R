normmix_marginal <- function(fit) {
  pattern <- cbind(fit$data$x,fit$data$y)
    den <-sapply( fit$allgens_List, densNormMix_atxy_sppmix, atxy = pattern)
    logden <- sapply(den,log)
    sumlogden <- sapply(logden,sum)
  marg <- mean(exp(sumlogden))
  return(marg)
}
