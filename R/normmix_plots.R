#' Plot normal mixture density.
#'
#' Plot the density surface of a normal mixture.
#'
#' @param intsurf An object of class \code{\link{intensity_surface}}.
#' @param win An object of class \code{\link[spatstat]{owin}}.
#' @param L Number of grid on x and y axes.
#' @param truncate Whether the density is truncated in the window (truncate)
#'  or not.
#' @param title1 The title (on top) for the density plot.
#' @param zlims The limits of z axis. The default does not has
#' additional limits on z axis.
#' @param poc The point coodinates for the title. The default is the top middle
#' in the 3-D plot.
#' @param ... Further arguments passed to \code{\link[rgl]{persp3d}}.
#'
#' @export
#' @examples
#' intsurf1 <- normmix(ps = c(.3, .7),
#'                     mus = list(c(0.2, 0.2), c(.8, .8)),
#'                     sigmas = list(.01*diag(2), .01*diag(2)),
#'                     lambda = 100,
#'                     win = square(1))
#' plot(intsurf1)
plot.intensity_surface <- function(intsurf, win = intsurf$window, L = 100,
                                   title1="Poisson Process Intensity",
                                   truncate = TRUE,
                                   zlims = c(0, 0),
                                   pos=c(0,0,0), grayscale = FALSE, ...) {
  xcoord <- seq(win$xrange[1], win$xrange[2], length.out = L)
  ycoord <- seq(win$yrange[1], win$yrange[2], length.out = L)

  xlims <- c(win$xrange)
  ylims <- c(win$yrange)

  z <- intsurf$intensity*dnormmix(intsurf, xlim = win$xrange, ylim = win$yrange,
                                  L = L, truncate = truncate)$v

  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                   "#7FFF7F", "yellow", "#FF7F00", "red",
                                   "#7F0000"))
  if (grayscale == TRUE) {
    col <- gray.colors(100)[findInterval(t(z), seq(min(z),
                                                   max(z), length = 100))]
  } else {
    col <- jet.colors(100)[findInterval(t(z), seq(min(z),
                                                  max(z), length = 100))]
  }

  if (zlims[1] == 0 && zlims[2] == 0) {
    zlims=c(0,max(z))
  }

  if (.Platform$OS.type == "unix") {
    rgl::layout3d(matrix(1:2, 1, 2), widths = c(5, 1))
  }
  rgl::open3d(windowRect = c(0, 45, 612, 657), zoom=1.2)
  U=rgl::par3d("userMatrix")
  rgl::par3d(userMatrix=
               rgl::rotate3d(U,pi/4,0,0,1))

  rgl::persp3d(x = xcoord, y = ycoord, z = t(z),
               color = col, xlab="x",ylab="y",zlab="",
               zlim=c(zlims[1]-0.01,zlims[2]),
               box = FALSE, axes = FALSE)
  rgl::axis3d('x')
  rgl::axis3d('y')
  rgl::axis3d('z-+', pos = c(xlims[1], ylims[2], 0))
  rgl::title3d(main = NULL)
  rgl::text3d(xlims[2], ylims[2], zlims[2], texts = title1)

  if (grayscale == TRUE) {
    rgl::bgplot3d(suppressWarnings(
      fields::image.plot(legend.only = TRUE,
                         #                         smallplot= c(.8,.82,0.05,.7),
                         zlim = zlims,
                         col = gray.colors(100))))
  } else {
    rgl::bgplot3d(suppressWarnings(
      fields::image.plot(legend.only = TRUE,
                         #                         smallplot= c(.8,.82,0.05,.7),
                         zlim = zlims,
                         col = jet.colors(100))))

  }
}


#' Plot spatial point pattern generated from mixture.
#'
#' Plot funciton for spatial point pattern generated from mixture. It's an
#' alternative to spatstat's plotting functions.
#'
#' @param pattern A point pattern of class sppmix or
#' \code{\link[spatstat]{ppp}}.
#' @param ... Parameters passed to \code{\link[ggplot2]{ggplot}}.
#' @inheritParams plot_contour
#' @export
#' @examples
#' intsurf1 <- normmix(ps = c(.3, .7),
#'                     mus = list(c(0.2, 0.2), c(.8, .8)),
#'                     sigmas = list(.01*diag(2), .01*diag(2)),
#'                     lambda = 100,
#'                     win = square(1))
#'
#' spp <- rsppmix(intsurf1)
#' plot(spp)
#' plot(spp, intsurf1)
#'
plot.sppmix <- function(pattern, mus, lambda,...) {
  n <- spatstat::npoints(pattern)

  p <- ggplot(as.data.frame(pattern), aes(x, y), ...) + geom_point() +
    labs(x = "X", y = "Y") +
    coord_fixed(xlim = spp$window$xrange, ylim = spp$window$yrange) +
    theme_light()

  if (missing(mus)) {
    p + ggtitle(paste("Spatial Point Pattern with", n, "points"))
  } else {
    stopifnot(is.intensity_surface(intsurf))
    mean_df <- data.frame(do.call(rbind, intsurf1$mus))
    names(mean_df) <- c("x", "y")
    p <- p + geom_point(data = mean_df, color = "red", size = 2.5) +
      ggtitle(paste("Spatial Point Pattern with", n, "points"))
  }
  p
}





#' 2D contour plots for normal mixture and intensity surface
#'
#' Create a contour plot for the intensity surface with or without realizations
#' on it.
#'
#' @inheritParams rsppmix
#' @param pattern optional spatial point pattern to add to the plot. In the c
#' class of \code{\link[spatstat]{ppp}}.
#' @param contour Logical flag indicating whether to plot countour only.
#' @param L number of grids on each coordinate.
#' @param ... Further parameters to passed to \code{to_int_surf()}.
#'
#' @import ggplot2
#' @examples
#' # plot normmix density
#' mix1 <- normmix(ps = c(.3, .7),
#'                 mus = list(c(0.2, 0.2), c(.8, .8)),
#'                 sigmas = list(.01*diag(2), .01*diag(2)))
#' plot(mix1)
#'
#' # plot intensity surface
#' intsurf1 <- normmix(ps = c(.3, .7),
#'                     mus = list(c(0.2, 0.2), c(.8, .8)),
#'                     sigmas = list(.01*diag(2), .01*diag(2)),
#'                     lambda = 100,
#'                     win = square(1))
#' plotmix_2d(intsurf1)
#'
#' pp1 <- rsppmix(intsurf = intsurf1)
#' plotmix_2d(intsurf1, pp1)  # with points
#'
#' @export
#' @rdname density_plots
plotmix_2d <- function(intsurf, pattern, contour = FALSE, truncate = TRUE,
                       L = 256, ...) {
  intsurf <- to_int_surf(intsurf, ...)
  win <- intsurf$window

  est_intensity <- dnormmix(intsurf, xlim = win$xrange, ylim = win$yrange,
                            L = L, truncate = truncate)

  p <- plot_density(as.data.frame(est_intensity), contour = contour) +
    labs(fill = "Intensity")
  if (!missing(pattern)) {
    p + geom_point(data=as.data.frame(pattern)) +
      ggtitle(bquote(paste("Normal Mixture Intensity Surface with ",
                           lambda == .(intsurf$intensity), ", ",
                           n == .(pattern$n))))
  } else {
    p + ggtitle(bquote(paste("Normal Mixture Intensity Surface with ",
                             lambda == .(intsurf$intensity))))
  }
}

#' @param mix object of class \code{normmix}.
#' @export
#' @rdname density_plots
plot.normmix <- function(mix, xlim, ylim, contour = FALSE,
                         truncate = FALSE, L = 256) {

  stopifnot(is.normmix(mix))

  if (missing(xlim) | missing(ylim)) {
    limits <- lapply(seq_len(mix$m), function(k) {
      c(mvtnorm::qmvnorm(.01, mean = mix$mus[[k]],
                         sigma = mix$sigmas[[k]])$quantile,
        mvtnorm::qmvnorm(.99, mean = mix$mus[[k]],
                         sigma = mix$sigmas[[k]])$quantile)
    })

    if (truncate) warning("Using truncation without specifying a window.")
    est_range <- range(unlist(limits))
    if (missing(xlim)) xlim <- est_range
    if (missing(ylim)) ylim <- est_range
  }

  est_density <- dnormmix(mix, xlim = xlim, ylim = ylim,
                          L = L, truncate = truncate)


  plot_density(as.data.frame(est_density), contour = contour) +
    labs(fill = "Density") +
    ggtitle("Normal Mixture Density Plot")
}


plot_density <- function(density_df, contour = FALSE) {
  p <- ggplot(density_df, aes(x, y)) + coord_fixed(expand = FALSE) +
    labs(x = "X", y = "Y")

  color <- c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
             "#FF7F00", "red", "#7F0000")

  if (!contour) {
    p + geom_raster(aes(fill = value)) +
      scale_fill_gradientn(colors = color) +
      guides(fill = guide_colorbar(nbin = 100, barheight = 15))
  } else {
    p + stat_contour(aes(z = value, color = ..level..)) +
      scale_color_gradientn(colors = color) +
      guides(color = guide_colorbar(nbin = 100, barheight = 15)) +
      theme_light()
  }
}








