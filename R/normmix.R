#' Mixture in 2d with normal components.
#'
#' Class of mixture in 2d with bivariate normal components.
#'
#' @param ps Vector of component probabilities.
#' @param mus A list where every element is a vector of length 2, defining the
#'  center of each component.
#' @param sigmas A list where every element is a 2 by 2 matrix, defining the
#'  covariance of each component.
#' @param lambda optional parameter of average intensity.
#' @param mix An object of class \code{\link{normmix}}
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
normmix <- function(ps, mus, sigmas, lambda) {
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

  if (!missing(lambda)) {
    if (lambda <= 0) stop("Intensity must be greater than 0.")
    RVAL <- c(RVAL, lambda = lambda)
    class(RVAL) <- c("normmix_intensity", "normmix")
  }

  RVAL
}

is.normmix <- function(mix) {
  # check class of object, if TRUE, assuming it's created by normmix() and
  # subject to all its constraints.
  ifelse(inherits(mix, "normmix"), TRUE, FALSE)
}

#' @rdname normmix
#' @export
print.normmix <- function(mix) {
  cat("Normal Mixture with", mix$m, "component(s).\n")
}


#' @rdname normmix
#' @export
print.normmix_intensity <- function(int_surf) {
  cat("Intensity surface with average intensity:", int_surf$lambda, "\n",
      "and distribution defined by: ")
  NextMethod(int_surf)
}


#' Generate mixture with normal components.
#'
#' Generate a mixture on a 2d window where the mean and variance of the
#' components are random. The number of component can either be fixed or random.
#'
#' @param m Number of component if \code{rand_m = FALSE}. When
#' \code{rand_m = TRUE}, number of component is random and this is the maximum
#' number of components.
#' @param sig0 Tunning parameter in generating random matrix from Wishart
#' distribution.
#' @param sigdf Degree of freedom in generating random matrix from Wishart
#' distribution.
#' @param win Object of class \code{spatstat::owin}. The mean vectors are
#' inside this window.
#' @param rand_m Whether the number of components are random or fixed (default).
#'
#' @return Object of class \code{normmix}.
#' @export
#' @examples
#' mix1 <- rnormmix(3, .01, 5, square(5))
#' mix2 <- rnormmix(8, .01, 10, square(1), rand_m = TRUE)
#'
rnormmix <- function(m, sig0, sigdf,
                     win, rand_m = FALSE) {
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
    sigmas[[k]] <- stats::rWishart(1, sigdf, sig0 * diag(2))[, , 1]
  }

  normmix(ps, mus, sigmas)
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

#' @rdname normmix
#' @export
summary.normmix_intensity <- function(normmix_intensity) {
  cat("Average number of points over domain:", normmix_intensity$lambda, "\n")
  NextMethod(normmix_intensity)
}




