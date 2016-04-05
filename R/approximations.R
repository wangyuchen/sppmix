#' Calculate density of normal mixture.
#'
#' When a \code{normmix} object is given, this function estimates density values
#' within the domain. When a \code{intensity_surface} is given, it multiplies
#' the estimated density values with intensity and returns the estimated
#' intensity.
#'
#' @param mix An object of class \code{normmix} or \code{intensity_surface}.
#' When an intensity surface is given, it estimates the intensity instead of
#' density within the domain.
#' @param xlim,ylim The range within which the density or intensity is
#' calculated.
#' @param L Length of grid on each axis. The density is calculated on a L * L
#' grid.
#' @param truncate Whether to truncate the density to be within the domain.
#'
#' @return An object of class \code{\link[spatstat]{im}}. This is a pixel image
#'  on a grid with values corresponding to the density at that location.
#'
#'
#' @examples
#' est_density <- dnormmix(demo_mix)
#' est_density <- dnormmix(demo_intsurf)
#' @export
dnormmix <- function(mix, xlim = c(0, 1), ylim = c(0, 1), L = 128,
                     truncate = TRUE) {

  if (!is.normmix(mix)) {
    stop("mix must be of class normmix or intensity surface")
  }

  if (is.intensity_surface(mix)) {
    xlim <- mix$window$xrange
    ylim <- mix$window$yrange
  }

  x <- seq(xlim[1], xlim[2], length.out = L)
  y <- seq(ylim[1], ylim[2], length.out = L)

  if (truncate) approx <- approx_normmix(mix, xlim, ylim)

  locs <- expand.grid(x, y)
  den <- matrix(NA_real_, nrow(locs), mix$m)
  for (k in 1:mix$m) {
    # every row of den is for a point
    # every col of den is for a component
    den[, k] <- mvtnorm::dmvnorm(locs, mix$mus[[k]], mix$sigmas[[k]])
    den[, k] <- den[, k] * mix$ps[k]

  }
  if (truncate) {
      den <- scale(den, center = FALSE, scale = approx)
  }

  est <- matrix(rowSums(den), nrow = length(y), ncol = length(x), byrow = T)

  if (is.intensity_surface(mix)) {
    spatstat::im(est * mix$intensity, x, y)
  } else {
    spatstat::im(est, x, y)
  }
}

