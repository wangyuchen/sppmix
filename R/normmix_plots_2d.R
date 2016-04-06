#' Plot spatial point pattern generated from mixture.
#'
#' Plot funciton for spatial point pattern generated from mixture. It's an
#' alternative to spatstat's plotting functions.
#'
#' @param x A point pattern of class sppmix or
#' \code{\link[spatstat]{ppp}}.
#' @param ... Additional parameters to \code{add_title} to be plot on title.
#' @examples
#' spp <- rsppmix(demo_intsurf)
#' plot(spp)
#' plot(spp, lambda = 100)
#'
#' @import ggplot2
#' @export
plot.sppmix <- function(x, ...) {
  n <- spatstat::npoints(x)
  win <- x$window

  p <- ggplot(as.data.frame(x), aes_string("x", "y")) + geom_point() +
    geom_rect(aes(xmin = win$xrange[1], xmax = win$xrange[2],
                  ymin = win$yrange[1], ymax = win$yrange[2]),
              color = "black", fill = NA) +
    labs(x = "X", y = "Y") +
    coord_fixed(xlim = win$xrange, ylim = win$yrange) +
    add_title("Spatial Point Pattern from Normal Mixture", n = x$n, ...)

  p
}


#' 2D density plots for normal mixture and intensity surface
#'
#' Create a contour plot for the intensity surface with or without realizations
#' on it.
#'
#' @inheritParams rsppmix
#' @param pattern Optional spatial point pattern to add to the plot. In the
#' class of \code{\link[spatstat]{ppp}}.
#' @param contour Logical flag indicating whether to plot countour only.
#' @param L Number of grids on each coordinate.
#'
#' @examples
#' # plot normmix density
#' plot(demo_mix)
#'
#' # plot intensity surface
#' plot(demo_intsurf)
#'
#' pp1 <- rsppmix(demo_intsurf)
#' plot(demo_intsurf, pattern = pp1)  # with points
#'
#' @export
#' @rdname density_plots
plot.intensity_surface <- function(x, ..., pattern, contour = FALSE,
                                   truncate = TRUE, L = 128) {
  x <- to_int_surf(x, ...)
  win <- x$window

  est_intensity <- dnormmix(x, xlim = win$xrange, ylim = win$yrange,
                            L = L, truncate = truncate)
  plot_df <- as.data.frame(est_intensity)
  names(plot_df) <- c("x", "y", "value")

  p <- plot_density(plot_df, contour = contour) +
    labs(fill = "Intensity")
  if (!missing(pattern)) {
    p + geom_point(data=as.data.frame(pattern)) +
      add_title("Normal Mixture Intensity Surface",
                lambda = x$intensity, m = x$m, n = pattern$n)
  } else {
    p + add_title("Normal Mixture Intensity Surface",
                  lambda = x$intensity, m = x$m)
  }
}

#' @param xlim,ylim Limits for the plot when \code{x} is a \code{normmix}. If
#' missing, quantile estimation is used to find a reasonable region.
#' @export
#' @rdname density_plots
plot.normmix <- function(x, xlim, ylim, contour = FALSE,
                         truncate = FALSE, L = 128, ...) {

  stopifnot(is.normmix(x))

  if (missing(xlim) | missing(ylim)) {
    limits <- lapply(seq_len(x$m), function(k) {
      c(mvtnorm::qmvnorm(.01, mean = x$mus[[k]],
                         sigma = x$sigmas[[k]])$quantile,
        mvtnorm::qmvnorm(.99, mean = x$mus[[k]],
                         sigma = x$sigmas[[k]])$quantile)
    })

    if (truncate) warning("Using truncation without specifying a window.")
    est_range <- range(unlist(limits))
    if (missing(xlim)) xlim <- est_range
    if (missing(ylim)) ylim <- est_range
  }

  est_density <- dnormmix(x, xlim = xlim, ylim = ylim,
                          L = L, truncate = truncate)
  plot_df <- as.data.frame(est_density)
  names(plot_df) <- c("x", "y", "value")


  plot_density(plot_df, contour = contour) +
    labs(fill = "Density") +
    add_title("Normal Mixture Density Plot", m = x$m)
}


plot_density <- function(density_df, contour = FALSE) {
  p <- ggplot(density_df, aes_string("x", "y")) + coord_fixed(expand = FALSE) +
    labs(x = "X", y = "Y")

  color <- c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
             "#FF7F00", "red", "#7F0000")

  if (!contour) {
    p + geom_raster(aes_string(fill = "value"), interpolate = TRUE) +
      scale_fill_gradientn(colors = color) +
      guides(fill = guide_colorbar(nbin = 100, barheight = 15))
  } else {
    p + stat_contour(aes_string(z = "value", color = "..level..")) +
      scale_color_gradientn(colors = color) +
      guides(color = guide_colorbar(nbin = 100, barheight = 15))
  }
}


add_title <- function(title, lambda = "", m = "", n = "", L = "") {
  if (!missing(lambda)) lambda <- bquote(paste(lambda == .(lambda)))
  if (!missing(m)) m <- bquote(paste(m == .(m), " components"))
  if (!missing(n)) n <- bquote(paste(n == .(n), " points"))
  if (!missing(L)) L <- bquote(paste(L == .(L), " iterations"))

  all_char <- list(lambda = lambda, m = m, n = n, L = L)
  non_empty_char <- all_char[nchar(all_char) > 0]

  cal <- do.call(function(...) substitute(list(...)), non_empty_char)

  ggtitle(substitute(atop(title, cal)))
}







