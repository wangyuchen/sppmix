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

#' Plot contour for intensity surface
#'
#' Create a contour plot for the intensity surface with or without realizations
#' on it.
#'
#' @param intsurf Object of class \code{\link{intensity_surface}}.
#' @param pattern Object of class \code{\link{ppp}}
#' @param win Object of class \code{\link{owin}}
#' @param L number of grids on each coordinate.  The default is L=100.
#' @param title title for the contour plot.
#' @param filled Logical flag indicating whether plot filled contour plot.
#' The default is TRUE.
#' @param truncate Logical flag indicating whether truncation is used for
#'  \code{pattern}. The default is TRUE.
#' @param ... Further arguments passed to \code{\link[rgl]{filled.contour}}.
#'
#' @examples
#' intsurf1 <- normmix(ps = c(.3, .7),
#'                     mus = list(c(0.2, 0.2), c(.8, .8)),
#'                     sigmas = list(.01*diag(2), .01*diag(2)),
#'                     lambda = 100,
#'                     win = square(1))
#' pp1 <- rsppmix(intsurf = intsurf1)
#'
#' # plot the theoretical intensity surface with a realization
#' plot_contour(intsurf1, pp1, spatstat::square(1))
#' @export
plot_contour <- function(intsurf, pattern, win = intsurf$window, L = 256,
                         title="Mixture intensity with", filled = TRUE, truncate = TRUE, ...) {
  if(!missing(pattern)) {n <- pattern$n}
  lambda <- intsurf$intensity
  xcoord <- seq(win$xrange[1], win$xrange[2], length.out = L)
  ycoord <- seq(win$yrange[1], win$yrange[2], length.out = L)

  surf <- dnormmix(intsurf, win = win, xcoord, ycoord, truncate = truncate)
  z <- intsurf$intensity*surf$v
  grid <- expand.grid(xcoord,ycoord)
  temp <- data.frame(grid$Var1,grid$Var2,as.vector(t(z)))
  names(temp) <- c("x","y","z")

  # assign colors to heights for each point
  jet.colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                              "#7FFF7F", "yellow", "#FF7F00", "red",
                                              "#7F0000"))
  #col <- jet.colors(1000)[findInterval(z$v, seq(min(z$v), max(z$v), length = 1000))]

  #plot(z, col = col)
  if (filled == TRUE){
    obj <- list(x = xcoord, y= ycoord, z= t(z))
    fields::plot.surface(obj, type = "I", col=jet.colors(150), main = "", ...)
    if (!missing(pattern)){
      points(pattern$x,pattern$y, pch=16, cex=0.5)
      title1 <- list(
        bquote(paste(lambda,"=",.(lambda),", n=",.(n))), title)
      mtext(do.call(expression, title1), side = 3, line = 0:1, at = 0.5)
    } else {
      title2 <- list(bquote(paste(lambda,"=",.(lambda))), title)
      mtext(do.call(expression, title2), side = 3, line = 0:1, at = 0.5)
    }
  } else {
    m <- ggplot2::ggplot(temp,aes(x, y, z = z)) +
      ggplot2::labs(x = "X", y = "Y",
                    title = paste("Contour of the intensity surface with",
                                  pattern$n, "Points")) +
      ggplot2::stat_contour(aes(colour = ..level..))
    if (!missing(pattern)) m <- m + geom_point(data=as.data.frame(pattern),aes(x, y, z=0))
    return(m)
  }
}

