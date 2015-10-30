#' @export
plot_avgsurf <- function(gens, LL = 30, burnin = 1000, xlims = c(0,1),
                                ylims = c(0,1)) {

  mix_of_postmeans <- MakeMixtureList(gens$allgens_List,burnin)
  mean_lambda <- mean(gens$genlamdas[burnin:L])
  zmax_genmeanmix <- mean_lambda * GetMixtureMaxz_sppmix(mix_of_postmeans,
                          LL,xlims,ylims)
  #find the highest z
  maxz_height <- max(zmax_genmeanmix)
  zlims <- c(0, 1.1*maxz_height)
  L  <-  dim(gens$genps)[1]
  gridvals  <-  GetGrid_sppmix(LL,xlims,ylims);
  ticsx <- gridvals[[1]];
  ticsy <- gridvals[[2]];
    ApproxAvgPostIntensityz <- ApproxAvgPostIntensity(
    gens$allgens_List,gens$genlamdas,LL,burnin,
    xlims,ylims)

  Plot3d_sppmix(xcoord = as.vector(ticsx),
                ycoord = as.vector(ticsy),
                zcoord = ApproxAvgPostIntensityz,
                title1 = paste("Average of",L-burnin,"posterior realizations
                             of the intensity function")
                ,xlims,ylims,zlims)


}
