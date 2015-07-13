#' Estimating Riemann integral over a 2d doamin.
#'
#' Estimating Riemann integral of a function over a 2d doamin. The domain is partitioned into small regions and the function is evaluated on the center of each region.
#'
#' @param func A vectorized function with two inputs. It should take two vector of the same length and return a numerical vector.
#' @param xlimits A vector with two values defining the range of the domain on x-axis.
#' @param ylimits A vector with two values defining the range of the domain on y-axis.
#' @param L A numerical value controls how many sub-regions are partitioned on each axis, default to 100.
#' @return The estimated value of Riemann integral of the function on that domain.
#' @examples
#' func <- function(x, y) exp(-.5 * x^2 - .5 * y^2) / (2 * pi)
#' est_riemann_int(func, c(-3, 3), c(-3, 3))
est_riemann_int <- function(func, xlimits, ylimits, L=100) {
  # lower leverl function for estimating Riemann integral

  if (length(xlimits) != 2 | length(ylimits) != 2) {
    stop("Limits must be vector of length 2.")
  }

  if (xlimits[1] >= xlimits[2] | ylimits[1] >= ylimits[2]) {
    stop("Lower limits must be smaller than upper limits")
  }

  x <- seq(from=xlimits[1], to=xlimits[2], length.out=L)
  y <- seq(from=ylimits[1], to=ylimits[2], length.out=L)

  lengths <- expand.grid(x_length=diff(x),
                         y_length=diff(y))
  centers <- expand.grid(x_center=(x[1:(L-1)] + x[2:L])/2,
                         y_center=(y[1:(L-1)] + y[2:L])/2)

  return(sum(lengths$x_length * lengths$y_length *
               func(centers$x_center, centers$y_center)))
}


#' Approximate density of a normal mixture over a 2d domain.
#'
#' Approximate the density of each component in a normal mixture within the domain using Riemann integral.
#'
#' @param mix An object of class normmix
#' @param win An object of class owin
#'
#' @return A numerical vector corresponding to the density of each component within thw window.
#'
#' @examples
#' mix1 <- normmix(ps=c(.5, .5),
#'                 mus=list(c(0, 0), c(1, 1)),
#'                 sigmas=list(.01*diag(2), .01*diag(2)))
#' approx_normmix(mix1, spatstat::square(1))
approx_normmix <- function(mix, win) {
  if (!is.normmix(mix)) {
    stop("mix must be an object of class normmix.")
  }

  if (!spatstat::is.owin(win)) {
    stop("win must be of class owin.")
  }

  approx <- numeric(mix$m)

  for (k in 1:mix$m) {
    func <- function(x, y) mvtnorm::dmvnorm(cbind(x, y),
                                            mix$mus[[k]], mix$sigmas[[k]])
    approx[k] <- est_riemann_int(func, win$xrange, win$yrange)
  }
  return(approx)
}







