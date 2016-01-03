#' Plot the average intensity surface
#'
#' Plot the average of the posterior realizations for the intensity surface
#'
#' @param fit An object contains all posterior realizations from
#' \code{\link{est_mix_damcmc}} or \code{\link{est_mix_bdmcmc}}.
#' @param win An object of class \code{\link[spatstat]{owin}}.
#' @param LL Number of grid on x and y axes.
#' @param burnin Length of burnin, default value is 1/10 of total number of
#' iteration.
#' @param zlims The limits of z axis. The default does not has
#' additional limits on z axis.
#' @author Athanasios Christou Micheas, Jiaxun Chen, Yuchen Wang
#' @export
#' @examples
#' # generate data
#' mix2 <- normmix(ps=c(.4, .6), mus=list(c(0.1, 0.1), c(0.8, 0.8)),
#' sigmas=list(.02*diag(2), .01*diag(2)))
#' pp2 <- rsppmix(100,mix2,square(1))
#' # Run Data augmentation MCMC and get posterior realizations
#' post=est_mix_damcmc(pp2,L = 5000,2,truncate = F)
#' # Plot the average of realized surfaces
#' plot_avgsurf(fit = post, win = square(1), LL = 30, burnin = 1000)
plot_avgsurf <- function(fit, win, LL = 30,
                         burnin = length(fit$allgens_List) / 10,
                         zlims = c(0, 0)) {

  # get limits
  xlims <- c(win$xrange)
  ylims <- c(win$yrange)
  L  <-  dim(fit$genps)[1]

  mix_of_postmeans <- MakeMixtureList(fit$allgens_List,burnin)
  mean_lambda <- mean(fit$genlamdas[burnin:L])
  zmax_genmeanmix <- mean_lambda * GetMixtureMaxz_sppmix(mix_of_postmeans,
                                                         LL,xlims,ylims)
  #find the highest z
  maxz_height <- max(zmax_genmeanmix)
  if (zlims[1] == 0 && zlims[2] == 0) {
    zlims <- c(0, 1.1*maxz_height)
  }
  gridvals  <-  GetGrid_sppmix(LL,xlims,ylims)
  xcoord <- as.vector(gridvals[[1]])
  ycoord <- as.vector(gridvals[[2]])
    zcoord <- ApproxAvgPostIntensity(
    fit$allgens_List, fit$genlamdas, LL, burnin,
    xlims, ylims)
  title1 = paste("Average of",L-burnin,
                 "posterior realizations of the intensity function")
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                   "#7FFF7F", "yellow", "#FF7F00", "red",
                                   "#7F0000"))
  col <- jet.colors(100)[findInterval(zcoord, seq(min(zcoord), max(zcoord),
                                                  length = 100))]

  #   if(.Platform$OS.type=='windows')
  #   {
  #     scr_width <- system("wmic desktopmonitor get screenwidth", intern=TRUE)
  #     scr_height <- system("wmic desktopmonitor get screenheight", intern=TRUE)
  #     height=as.numeric(scr_height[length(scr_height)-1])
  #     width=as.numeric(scr_width[length(scr_width)-1])
  #   }

  #  rgl::layout3d(matrix(1:2, 1, 2), widths = c(5, 1))
  rgl::open3d(windowRect = c(0, 45, 612, 657), zoom=1.2)

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
  rgl::bgplot3d(suppressWarnings(
    fields::image.plot(legend.only = TRUE,
                       smallplot= c(.8,.82,0.05,.7),
                       zlim = zlims,
                       col = jet.colors(100))))
}
