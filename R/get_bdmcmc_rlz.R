#' @export
get_bdmcmc_rlz <- function(fit, comp) {

  L <- length(fit$allgens)

  res <- list()
  ind <- fit$numcomp == comp

  res$allgens_List <- fit$allgens_List[ind]
  res$genps <- fit$genps[ind, seq_len(comp), drop = FALSE]
  res$genmus <- fit$genmus[seq_len(comp), , ind, drop = FALSE]

  res$gensigmas <- fit$gensigmas[ind, seq_len(comp), drop = FALSE]

  res$genlamdas <- fit$genlamdas[ind, 1, drop = FALSE]
  res$genzs <- fit$genzs[ind, , drop = FALSE]
  res$data <- fit$data

  class(res) <- "damcmc_res"
  return(res)
}
