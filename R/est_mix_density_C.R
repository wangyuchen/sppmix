#' @export
DAMCMC2d <- function(data,xlims, ylims, m, L, burnin, LL, trunc){
  .Call("sppmix_DAMCMC2d_sppmix",PACKAGE = "sppmix",data, xlims, ylims, m, L,
        burnin, LL, trunc)
}
