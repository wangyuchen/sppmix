est_riemann_int <- function(func, xlimits, ylimits, L=100) {
#   func <- function(x, y) exp(-.5 * x^2 - .5 * y^2) / (2 * pi)

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

