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
                                   title1="Poisson Process Intensity", truncate = TRUE,
                                   zlims = c(0, 0),
                                   pos=c(0,0,0), grayscale = FALSE, ...) {
  xcoord <- seq(win$xrange[1], win$xrange[2], length.out = L)
  ycoord <- seq(win$yrange[1], win$yrange[2], length.out = L)

  xlims <- c(win$xrange)
  ylims <- c(win$yrange)

  z <- intsurf$intensity*dnormmix(intsurf, win = win, xcoord, ycoord, truncate = truncate)$v

  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                   "#7FFF7F", "yellow", "#FF7F00", "red",
                                   "#7F0000"))
  if (grayscale == TRUE) {
    col <- gray.colors(100)[findInterval(t(z), seq(min(z), max(z), length = 100))]
  } else {
    col <- jet.colors(100)[findInterval(t(z), seq(min(z), max(z), length = 100))]
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
#' @param pattern A point pattern of class sppmix or \code{\link[spatstat]{ppp}}.
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
plot.sppmix <- function(pattern, intsurf, lambda,...) {
  n <- pattern$n
  if(!missing(intsurf)){
    if(is.intensity_surface(intsurf)) {
      lambda <- intsurf$intensity
      m <- intsurf$m
      title <- list(
        bquote(paste(lambda,"=",.(intsurf$intensity),",",.(n)," points, ", .(m),
                     " component(s)")), "Spatial Point Pattern with")
    }
  } else {
    title <- list(paste("Spatial Point Pattern with", n, "points"), NULL)
  }

  plot(pattern$x,pattern$y, xlab = "X", ylab = "Y", pch=16,main = "",...)
  mtext(do.call(expression, title),side=3,line=0:1)
  if (!missing(intsurf)) {
    center <- as.data.frame(matrix(unlist(intsurf$mus),ncol=2, byrow = T))
    colnames(center) <- c("x", "y")
    points(center$x, center$y, col = "red", pch=19)
  }
}

#' 2D plot for intensity surface
#'
#' Create a contour plot for the intensity surface with or without realizations
#' on it.
#'
#' @inheritParams rsppmix
#' @param pattern optional spatial point pattern to add to the plot. In the c
#' class of \code{\link[spatstat]{ppp}}.
#' @param filled Logical flag indicating whether plot filled contour plot.
#' The default is TRUE.
#' @param L number of grids on each coordinate.
#' @param density whether to plot mixture density instead of intensity.
#' @param xlim,ylim omitted except when density = TRUE
#'
#' @import ggplot2
#' @examples
#' intsurf1 <- normmix(ps = c(.3, .7),
#'                     mus = list(c(0.2, 0.2), c(.8, .8)),
#'                     sigmas = list(.01*diag(2), .01*diag(2)),
#'                     lambda = 100,
#'                     win = square(1))
#' pp1 <- rsppmix(intsurf = intsurf1)
#'
#' # plot the theoretical intensity surface with a realization
#' plotmix_2d(intsurf1, pp1)
#' @export
plotmix_2d <- function(intsurf, pattern, filled = TRUE, truncate = TRUE,
                       L = 256, density = FALSE, xlim, ylim, ...) {


  intsurf <- to_int_surf(intsurf, ..., return_normmix = density)

  win <- if (density) list(xrange = xlim, yrange = ylim) else intsurf$window

  est_den <- dnormmix(intsurf, xlim = win$xrange, ylim = win$yrange,
                      L = L, truncate = truncate)

  p <- ggplot(as.data.frame(est_den), aes(x, y)) +
    coord_fixed(xlim = win$xrange, ylim = win$yrange, expand = FALSE) +
    labs(x = "X", y = "Y")

  if (filled) {
    color <- c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
               "#FF7F00", "red", "#7F0000")

    p <- p + geom_raster(aes(fill = value)) +
      scale_fill_gradientn(colors = color) +
      guides(fill = guide_colorbar(nbin = 100, barheight = 15))

    common_title <- "Normal Mixture"

    if (density) {
      p <- p + ggtitle(paste(common_title, "Density Plot"))
    } else {
      p <- p + ggtitle(bquote(paste(.(common_title), "Intensity Surface with ",
                                    lambda == .(intsurf$intensity))))
    }

  } else {
    common_title <- "Contour of Normal Mixture"
    p <- p + stat_contour(aes(z = value))

    if (density) {
      p <- p + ggtitle(paste(common_title, "Density"))
    } else {
      p <- p + ggtitle(bquote(paste(.(common_title), "Intensity Surface with ",
                                    lambda == .(intsurf$intensity))))
    }
  }

  if (!missing(pattern)) {
    p <- p + geom_point(data=as.data.frame(pattern))
  }

  p
}

