plot_density_3d <- function(x, y, z, main, ...) {
  sppmix::demo_intsurf
  win <- demo_intsurf$window
  est_intensity <- dnormmix(demo_intsurf, xlim = win$xrange, ylim = win$yrange,
                            L = 64)
  est_intensity

  x <- est_intensity$xcol
  y <- est_intensity$yrow
  z <- as.matrix(est_intensity)

  # above should be passed in

  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                   "#7FFF7F", "yellow", "#FF7F00", "red",
                                   "#7F0000"))

  col <- jet.colors(100)[findInterval(z, seq(min(z), max(z), length = 100))]


  rgl::open3d(useNULL = TRUE)  # plot at backend when calling any *3d functions
  rgl::persp3d(x, y, z, col = col)

  s <- rgl::scene3d()
  class(s) <- c("density_3d", class(s))
  s
}

print.density_3d <- function(scene, new_window = FALSE, ...) {
  rglwidget::rglwidget(scene, ...)
}
