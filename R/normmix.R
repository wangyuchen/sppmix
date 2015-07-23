#' Mixture in 2d with normal components.
#'
#' Class of mixture in 2d with bivariate normal components.
#'
#' @param ps Vector of component probabilities.
#' @param mus A list where every element is a vector of length 2, defining the
#'  center of each component.
#' @param sigmas A list where every element is a 2 by 2 matrix, defining the
#'  covariance of each component.
#' @param mix An objecto of class \code{\link{normmix}}
#'
#' @return An object of class "normmix" containing the following components:
#'  \item{m}{Number of components.}
#'  \item{ps}{Vector of component probabilities.}
#'  \item{mus}{List of mean vectors of components.}
#'  \item{sigmas}{List of covariance matrix of components.}
#' @seealso \code{\link{rnormmix}} for generating random mixture.
#' @export
#' @examples
#' mix1 <- normmix(ps=c(.3, .7), mus=list(c(0, 0), c(1, 1)),
#'                 sigmas=list(.01*diag(2), .01*diag(2)))
#' mix1
#' summary(mix1)
#' if (require(spatstat)) {
#'   plot(mix1, square(1))
#' }
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

#' @rdname normmix
#' @export
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
#' @export
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

#' Calculate density of locations under a mixture.
#'
#' Calculate the density values of a given mixture with normal components.
#'
#' @inheritParams plot.normmix
#' @param x,y Locations where the density function will be evaluated. It will
#'  automatically expand it into a grid. If missing, a regular grid over the
#'  window will be given.
#'
#' @return An object of class \code{\link[spatstat]{im}}. This is a pixel image
#'  on a grid with values corresponding to the density at that location.
#'
#' @seealso \code{\link[spatstat]{summary.im}} and
#'  \code{\link[spatstat]{plot.im}} for manipulating with pixel image object.
#' @export
#' @examples
#' if (require(spatstat)) {
#'   mix1 <- rnormmix(8, sig0 = .01, 10, square(2))
#'   den <- dnormmix(mix1, square(2))
#'  plot(den)
#' }
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

#' Calculate intensity of locations under a mixture.
#'
#' Calculate intensity values of a given mixture at given locations.
#'
#' @param lambda Mean intensity value
#' @inheritParams dnormmix
#'
#' @return A pixel image of class \code{\link[spatstat]{im}} with the values
#'  corresponding to the intensity at given locations.
#' @export
#' @examples
#' if (require(spatstat)) {
#'   mix1 <- rnormmix(8, sig0 = .01, 10, square(2))
#'   int <- inormmix(100, mix1, square(2))
#'  plot(int)
#' }
inormmix <- function(lambda, mix, win, x, y, L = 128, truncate = TRUE) {
  if (missing(x)) x <- seq(win$xrange[1], win$xrange[2], length.out = L)
  if (missing(y)) y <- seq(win$yrange[1], win$yrange[2], length.out = L)

  intensity <- lambda * dnormmix(mix, win, x, y,
                                 truncate = truncate)$v

  RVAL <- spatstat::im(matrix(intensity, nrow = length(y),
                              ncol = length(x), byrow = T), x, y)
  return(RVAL)
}

#' @rdname normmix
#' @export
summary.normmix <- function(mix) {
  cat(paste("Normal Mixture with", mix$m, "components:\n"))
  for (i in 1:mix$m) {
    cat(paste("Component", i, "is centered at",
              "(", round(mix$mus[[i]][1], 2), ",",
              round(mix$mus[[i]][2], 2), ")",
              "with probability", round(mix$ps[[i]], 4), "\n"))
  }
}


as.normmix <- function(ps, mus, sigmas) {
  # coerce mus and sigmas in DAMCMC format into normmix
  mu <- list()
  sigma <- list()
  for (i in 1:nrow(mus)) {
    mu[[i]] <- as.numeric(mus[i, ])
    sigma[[i]] <- sigmas[, , i]
  }
  return(normmix(ps, mu, sigma))
}

