#' Plot spatial point pattern generated from mixture.
#'
#' Plot funciton for spatial point pattern generated from mixture. It's an
#' alternative to spatstat's plotting functions.
#'
#' @param pattern A point pattern of class sppmix or
#' \code{\link[spatstat]{ppp}}.
#' @param mus An optioanl list of theoretical mean for each component.
#' @param ... Additional parameters to \code{add_title} to be plot on title.
#' @examples
#' spp <- rsppmix(demo_intsurf)
#' plot(spp)
#' plot(spp, mus = intsurf1$mus)
#' plot(spp, mus = intsurf1$mus, lambda = 100)
#'
#' @import ggplot2
#' @export
plot.sppmix <- function(pattern, mus, ...) {
  n <- spatstat::npoints(pattern)

  p <- ggplot(as.data.frame(pattern), aes(x, y)) + geom_point() +
    labs(x = "X", y = "Y") +
    coord_fixed(xlim = pattern$window$xrange, ylim = pattern$window$yrange) +
    theme_light() +
    add_title("Spatial Point Pattern from Normal Mixture", n = pattern$n, ...)


  if (!missing(mus)) {
    mean_df <- data.frame(do.call(rbind, mus))
    names(mean_df) <- c("x", "y")
    p <- p + geom_point(data = mean_df, color = "red", size = 2.5) +
      add_title("Spatial Point Pattern from Normal Mixture", n = pattern$n,
                m = nrow(mean_df), ...)
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
#'
#' @examples
#' # plot normmix density
#' plot(demo_mix)
#'
#' # plot intensity surface
#' plotmix_2d(demo_intsurf)
#'
#' pp1 <- rsppmix(intsurf = demo_intsurf)
#' plotmix_2d(demo_intsurf, pp1)  # with points
#'
#' @export
#' @rdname density_plots
plot.intensity_surface <- function(intsurf, pattern, contour = FALSE,
                                   truncate = TRUE, L = 256, ...) {
  intsurf <- to_int_surf(intsurf, ...)
  win <- intsurf$window

  est_intensity <- dnormmix(intsurf, xlim = win$xrange, ylim = win$yrange,
                            L = L, truncate = truncate)

  p <- plot_density(as.data.frame(est_intensity), contour = contour) +
    labs(fill = "Intensity")
  if (!missing(pattern)) {
    p + geom_point(data=as.data.frame(pattern)) +
      add_title("Normal Mixture Intensity Surface",
                lambda = intsurf$intensity, m = intsurf$m, n = pattern$n)
  } else {
    p + add_title("Normal Mixture Intensity Surface",
                  lambda = intsurf$intensity, m = intsurf$m)
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
    add_title("Normal Mixture Density Plot", m = mix$m)
}


plot_density <- function(density_df, contour = FALSE) {
  p <- ggplot(density_df, aes(x, y)) + coord_fixed(expand = FALSE) +
    labs(x = "X", y = "Y")

  color <- c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
             "#FF7F00", "red", "#7F0000")

  if (!contour) {
    p + geom_raster(aes(fill = value), interpolate = TRUE) +
      scale_fill_gradientn(colors = color) +
      guides(fill = guide_colorbar(nbin = 100, barheight = 15))
  } else {
    p + stat_contour(aes(z = value, color = ..level..)) +
      scale_color_gradientn(colors = color) +
      guides(color = guide_colorbar(nbin = 100, barheight = 15)) +
      theme_light()
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







