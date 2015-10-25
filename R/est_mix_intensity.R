#' @export
est_mix_damcmc <- function(pp, m, truncate = FALSE,
                           L = 5000, burnin = 500, LL = 100) {

  fit <- DAMCMC2d_sppmix(data = cbind(pp$x, pp$y),
                         xlims = Window(pp)$xrange, ylims = Window(pp)$yrange,
                         m = m, truncate = truncate,
                         L = L, burnin = burnin, LL = LL)
  class(fit) <- "damcmc_res"
  return(invisible(fit))
}



#' @export
est_mix_bdmcmc <- function(pp, max_comp, truncate = FALSE,
                           lambda, lambdab, hyper,
                           L = 5000, burnin = 500, LL = 100) {

  fit <- BDMCMC2d_sppmix(max_comp, data = cbind(pp$x, pp$y),
                         xlims = Window(pp)$xrange, ylims = Window(pp)$yrange,
                         truncate = truncate,
                         lamda = lambda, lamdab = lambdab, hyper = hyper,
                         L = L, burnin = burnin, LL = LL)
  class(fit) <- "bdmcmc_res"
  return(invisible(fit))
}
