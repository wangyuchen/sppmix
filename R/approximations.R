#' Approximate density of a normal mixture over a 2d domain.
#'
#' Approximate the density of each component in a normal mixture within the
#' domain using multivariate normal density function.
#'
#' @param mix An object of class \code{\link{normmix}}
#' @param win An object of class \code{\link[spatstat]{owin}}
#'
#' @return A numerical vector corresponding to the density of each component
#'  within the window.
#' @export
#' @examples
#' if (require("spatstat")) {
#'   mix1 <- normmix(ps=c(.3, .7),
#'                   mus=list(c(0, 0), c(1, 1)),
#'                  sigmas=list(.01*diag(2), .01*diag(2)))
#'   win1 <- spatstat::square(1)
#'   plot(dnormmix(mix1, win1))
#'   approx_normmix(mix1, win1)
#' }
approx_normmix <- function(mix, win) {
  if (!is.normmix(mix)) {
    stop("mix must be an object of class normmix.")
  }

  if (!spatstat::is.owin(win)) {
    stop("win must be of class owin.")
  }

  approx <- numeric(mix$m)

  for (k in 1:mix$m) {
    approx[k] <- mvtnorm::pmvnorm(lower = c(win$xrange[1], win$yrange[1]),
                                  upper = c(win$xrange[2], win$yrange[2]),
                                  mean = mix$mus[[k]],
                                  sigma = mix$sigmas[[k]])
  }

  return(approx)
}

#' Calculate density of locations under a mixture.
#'
#' Calculate the density values of a given mixture with normal components.
#'
#' @inheritParams approx_normmix
#'
#' @return An object of class \code{\link[spatstat]{im}}. This is a pixel image
#'  on a grid with values corresponding to the density at that location.
#'
#' @seealso \code{\link[spatstat]{summary.im}} and
#'  \code{\link[spatstat]{plot.im}} for manipulating with pixel image object.
#' @export
dnormmix <- function(mix, xlim = c(0, 1), ylim = c(0, 1), L = 128,
                     truncate = TRUE) {

  if (!is.normmix(mix)) {
    stop("mix must be of class normmix or intensity surface")
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
    if (truncate) {
      den[, k] <- den[, k] / approx[k]
    }
  }

  est <- matrix(rowSums(den), nrow = length(y), ncol = length(x), byrow = T)

  if (is.intensity_surface(mix)) {
    spatstat::im(est * mix$intensity, x, y)
  } else {
    spatstat::im(est, x, y)
  }
}

