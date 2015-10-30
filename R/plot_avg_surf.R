#' Plot the average intensity surface
#'
#' Plot the average of the posterior realizations for the intensity surface
#'
#' @param gens An object contains all posterior realizations from
#' \code{\link{est_mix_damcmc}} or \code{\link{est_mix_bdmcmc}}.
#' @param win An object of class \code{\link[spatstat]{owin}}.
#' @param LL Number of grid on x and y axes.
#' @param burnin Length of burnin.
#'
#' @export
#' @examples
#' # generate data
#' mix2 <- normmix(ps=c(.4, .6), mus=list(c(0.1, 0.1), c(0.8, 0.8)),
#' sigmas=list(.02*diag(2), .01*diag(2)))
#' pp2 <- rsppmix(100,mix2,square(1))
#' # Run Data augmentation MCMC and get posterior realizations
#' post=est_mix_damcmc(pp2,L = 5000,2,truncate = F)
#' # Plot the average of realized surfaces
#' plot_avgsurf(gens = post, win = square(1), LL = 30, burnin = 1000)

plot_avgsurf <- function(gens, win, LL = 30, burnin = 1000) {

  # get limits
  xlims <- c(win$xrange)
  ylims <- c(win$yrange)

  mix_of_postmeans <- MakeMixtureList(gens$allgens_List,burnin)
  mean_lambda <- mean(gens$genlamdas[burnin:L])
  zmax_genmeanmix <- mean_lambda * GetMixtureMaxz_sppmix(mix_of_postmeans,
                          LL,xlims,ylims)
  #find the highest z
  maxz_height <- max(zmax_genmeanmix)
  zlims <- c(0, 1.1*maxz_height)
  L  <-  dim(gens$genps)[1]
  gridvals  <-  GetGrid_sppmix(LL,xlims,ylims)
  xcoord <- as.vector(gridvals[[1]])
  ycoord <- as.vector(gridvals[[2]])
    zcoord <- ApproxAvgPostIntensity(
    gens$allgens_List, gens$genlamdas, LL, burnin,
    xlims, ylims)
  title1 = paste("Average of",L-burnin,
                 "posterior realizations of the intensity function")
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                   "#7FFF7F", "yellow", "#FF7F00", "red",
                                   "#7F0000"))
  col <- jet.colors(100)[findInterval(zcoord, seq(min(zcoord), max(zcoord), length = 100))]

  if(is.null(zlims))
    zlims=c(0,max(zcoord))
  height=300;
  width=500;
#   if(.Platform$OS.type=='windows')
#   {
#     scr_width <- system("wmic desktopmonitor get screenwidth", intern=TRUE)
#     scr_height <- system("wmic desktopmonitor get screenheight", intern=TRUE)
#     height=as.numeric(scr_height[length(scr_height)-1])
#     width=as.numeric(scr_width[length(scr_width)-1])
#   }

  rgl::open3d(windowRect=c(width/5,
                           height/7,
                           4*width/5,
                           6*height/7),
              zoom=1.2)
  U=rgl::par3d("userMatrix")
  rgl::par3d(userMatrix=
               rgl::rotate3d(U,pi/4,0,0,1))
  zmax=max(zcoord)
  Rangez=zmax-min(zcoord);
  rgl::persp3d(x = xcoord, y = ycoord, z = zcoord,
               color = col, xlab="x",ylab="y",zlab="",
               zlim=c(zlims[1]-0.01,zlims[2]),
               box = FALSE, axes = FALSE)
  rgl::axis3d('x')
  rgl::axis3d('y')
  rgl::axis3d('z-+', pos = c(xlims[1], ylims[2], 0))
  rgl::title3d(main=NULL)
  rgl::text3d(xlims[2],ylims[2],
              zlims[2]
              ,texts= title1)
}
