#' 3D density plots for normal mixture and intensity surface.
#'
#' Plot the density surface of a normal mixture in an interactive 3D plot.
#'
#' @inheritParams plot.intensity_surface
#' @param main Title for the density plot.
#'
#' @export
#' @examples
#' s <- plot_mix_3d(demo_intsurf)
#' s
#'
plot_mix_3d <- function(x, ..., truncate = TRUE, L = 128,
                        main = "Normal Mixture 3D Surface Plot") {
  x <- to_int_surf(x, ...)
  win <- x$window

  est_intensity <- dnormmix(x, xlim = win$xrange, ylim = win$yrange,
                            L = L, truncate = truncate)
  x <- est_intensity$xcol
  y <- est_intensity$yrow
  z <- as.matrix(est_intensity)

  plot_density_3d(x, y, z, main = main, ...)
}


plot_density_3d <- function(x, y, z, main, ...) {
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                   "#7FFF7F", "yellow", "#FF7F00", "red",
                                   "#7F0000"))(100)
  zcol <- cut(z, 100)

  # plot at backend when calling any *3d functions
  rgl::open3d(useNULL = TRUE,
              userMatrix = structure(c(0.882, 0.013, -0.472, 0, -0.472,
                                       0.054, -0.88, 0, 0.014, 0.998, 0.054,
                                       0, 0, 0, 0, 1), .Dim = c(4L, 4L)))

  rgl::layout3d(matrix(1:2, 2), heights = c(1, 10))
  rgl::text3d(max(x), max(y), max(z), texts = main)
  rgl::next3d()

  rgl::persp3d(x, y, z, col = jet.colors[zcol],
               xlab = "X", ylab = "Y", zlab = "", axes = FALSE)
  rgl::axes3d(c('x', 'y', 'z-+'))

  s <- rgl::scene3d()
  class(s) <- c("density_3d", class(s))
  s
}

#' @export
print.density_3d <- function(x, ...) {
  print(rglwidget::rglwidget(x, ...))
}
