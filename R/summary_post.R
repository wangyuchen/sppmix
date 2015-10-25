#' @export
summary.damcmc_res <- function(gens, dgt = 2) {

  m = dim(gens$genmus)[1]
  for (i in 1:m)
  {
    #true value and credible sets
    poststats = GetStats_sppmix(gens$genps[,i],alpha=0.05)
    poststats <- sapply(poststats, format, digits = dgt)

    cat("\n----------------Component ",i,"------------\n")
    #cat(paste("Probability: true =",truemix[[i]]$p))
    cat(paste("\nProbability: posterior mean =",
              poststats$Mean,"\n"))
    cat(paste(poststats$CredibleSetConfidence,
              "% Credible Set:\n[",poststats$CredibleSet[1],
              ",",poststats$CredibleSet[2],"]" ))

    poststats = GetStats_sppmix(gens$genmus[i,1,],alpha=0.05)
    poststats <- sapply(poststats, format, digits = dgt)

    #cat(paste("\nMean vector, x-coord: true =",truemix[[i]]$mu[1]))
    cat(paste("\nMean vector, x-coord: post mean =",
              poststats$Mean,"\n"))
    cat(paste(poststats$CredibleSetConfidence,
              "% Credible Set:\n[",poststats$CredibleSet[1],
              ",",poststats$CredibleSet[2],"]" ))

    poststats = GetStats_sppmix(gens$genmus[i,2,],alpha=0.05)
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
