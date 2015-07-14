normmix <- function(ps, mus, sigmas) {
  if (length(ps) == 0) {
    stop("Mixture with 0 components")
  }

  if (abs(sum(ps) - 1) > .0000001) {
    stop("Component probabilities must sum to 1.")
  }

  if (length(ps) != length(mus) | length(ps) != length(sigmas)) {
    stop("Number of components mismatch.")
  }

  RVAL <- list(m = length(ps), ps = ps, mus = mus, sigmas = sigmas)
  class(RVAL) <- "normmix"
  return(RVAL)
}

is.normmix <- function(mix) {
  # check class of object, if TRUE, assuming it's created by normmix() and
  # subject to all its constraints.
  return(ifelse(class(mix) == "normmix", TRUE, FALSE))
}

print.normmix <- function(mix) {
  print(paste("Normal Mixture with", mix$m, "component"))
}

#' Generate mixture with normal components.
#'
#' Generate a mixture on a 2d window where the mean and variance of the components are random. The number of component can either be fixed or random.
#'
#' @param m Number of component if \code{rand_m = FALSE}. When \code{rand_m = TRUE}, number of component is random and this is the maximum number of components.
#' @param sig0 Tunning parameter in generating random matrix from Wishart distribution.
#' @param sigdf Degree of freedom in generating random matrix from Wishart distribution.
#' @param win Object of class \code{spatstat::owin}. The mean vectors are inside this window.
#' @param rand_m Whether the number of components are random or fixed (default).
#'
#' @return Object of class \code{normmix}.
#'
#' @examples
#' mix1 <- rnormmix(3, .01, 5, square(5))
#' mix2 <- rnormmix(8, .01, 10, square(1), rand_m = TRUE)
#'
rnormmix <- function(m, sig0, sigdf,
                     win,
                     rand_m=FALSE) {
  if (!spatstat::is.owin(win)) {
    stop("win must be of class owin.")
  }
  if (rand_m == TRUE) {
    # number of components is random
    m <- sample(1:m, 1)
  }

  ps <- rdirichlet(1, rep(1, m))
  mus <- list()
  sigmas <- list()

  for (k in 1:m) {
    mus[[k]] <- c(runif(1, win$xrange[1], win$xrange[2]),
                  runif(1, win$yrange[1], win$yrange[2]))
    sigmas[[k]] <- rWishart(1, sigdf, sig0 * diag(2))[, , 1]
  }

  return(normmix(ps, mus, sigmas))
}

dnormmix <- function(mix, win, x, y, L = 128, truncate = TRUE) {
  if (missing(x)) x <- seq(win$xrange[1], win$xrange[2], length.out = L)
  if (missing(y)) y <- seq(win$yrange[1], win$yrange[2], length.out = L)

  locs <- expand.grid(x, y)
  approx <- approx_normmix(mix, win)
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

  RVAL <- spatstat::im(matrix(rowSums(den), nrow = length(y),
                              ncol = length(x), byrow = T), x, y)
  return(RVAL)
}



inormmix <- function(lambda, mix, win, x, y, L = 128, truncate = TRUE) {
  if (missing(x)) x <- seq(win$xrange[1], win$xrange[2], length.out = L)
  if (missing(y)) y <- seq(win$yrange[1], win$yrange[2], length.out = L)

  intensity <- lambda * dnormmix(mix, win, x, y,
                                 truncate = truncate)$v

  RVAL <- spatstat::im(matrix(intensity, nrow = length(y),
                              ncol = length(x), byrow = T), x, y)
  return(RVAL)
}

summary.normmix <- function(mix) {
  cat(paste("Normal Mixture with", mix$m, "components:\n"))
  for (i in 1:mix$m) {
    cat(paste("Component", i, "is centered at",
              "(", round(mix$mus[[i]][1], 2), ",",
              round(mix$mus[[i]][2], 2), ")",
              "with probability", round(mix$ps[[i]], 4), "\n"))
  }
}

plot.normmix <- function(mix, win, L = 20, truncate=TRUE) {
  rgl::open3d()
  xcoord <- seq(win$xrange[1], win$xrange[2], length.out = L)
  ycoord <- seq(win$yrange[1], win$yrange[2], length.out = L)

  z <- dnormmix(mix, win = win, xcoord, ycoord, truncate = truncate)$v
  rgl::surface3d(xcoord, ycoord, z,
                 color="#FF2222", alpha=0.5,
                 rgl::axes3d(edges ="bbox", box=F, labels = T))

}
