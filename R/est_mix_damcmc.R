#' @export
est_mix_damcmc <- function(pp, m, truncate = TRUE,
                           L = 50000, burnin = 500, LL = 100) {

  fit <- .Call('sppmix_DAMCMC2d_sppmix', PACKAGE = 'sppmix',
               data = cbind(pp$x, pp$y),
               xlims = Window(pp)$xrange, ylims = Window(pp)$yrange,
               m = m, truncate = truncate,
               L = L, burnin = burnin, LL = LL)
  class(fit) <- "damcmc_res"
  return(fit)
}
