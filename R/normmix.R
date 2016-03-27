#' Mixture in 2d with normal components.
#'
#' Class of mixture in 2d with bivariate normal components.
#'
#' @param ps Vector of component probabilities.
#' @param mus A list where every element is a vector of length 2, defining the
#'  center of each component.
#' @param sigmas A list where every element is a 2 by 2 covariance matrix,
#' defining the covariance for each component.
#' @param lambda optional parameter of theoretical average number of points.
#' If set, the returned object will be an intensity surface.
#' @param win optional parameter of domain. Must be set with intensity value to
#' create an intensity surface.
#'
#' @return An object of class "normmix" containing the following components:
#'  \item{m}{Number of components.}
#'  \item{ps}{Vector of component probabilities.}
#'  \item{mus}{List of mean vectors of components.}
#'  \item{sigmas}{List of covariance matrix of components.}
#'  \item{intensity}{optional intensity value if lambda is provided when calling
#'  \code{normmix}.}
#'  \item{window}{optional window object of class \code{\link{spatstat::owin}}.}
#'
#' @seealso \code{\link{rnormmix}} for generating random mixture.
#' @export
#' @examples
#' mix1 <- normmix(ps = c(.3, .7),
#'                 mus = list(c(0.2, 0.2), c(.8, .8)),
#'                 sigmas = list(.01*diag(2), .01*diag(2)))
#' mix1
#'
#'
#' intsurf1 <- normmix(ps = c(.3, .7),
#'                     mus = list(c(0.2, 0.2), c(.8, .8)),
#'                     sigmas = list(.01*diag(2), .01*diag(2)),
#'                     lambda = 100,
#'                     win = square(1))
#' intsurf1
#'
normmix <- function(ps, mus, sigmas, lambda = NULL, win = NULL,
                    estimated = FALSE) {
  if (abs(sum(ps) - 1) > .001) {
    stop("Component probabilities must sum to 1.")
  }

  if (length(ps) != length(mus) | length(ps) != length(sigmas)) {
    stop("Number of components mismatch.")
  }

  if (!is.null(lambda) & spatstat::is.owin(win)) {
    # generating intensity surface
    if (lambda <= 0) stop("Intensity must be greater than 0.")
    RVAL <- list(m = length(ps), ps = ps, mus = mus, sigmas = sigmas,
                 intensity = lambda, window = win, estimated = estimated)
    class(RVAL) <- c("intensity_surface", "normmix")
  } else {
    RVAL <- list(m = length(ps), ps = ps, mus = mus, sigmas = sigmas,
                 estimated = estimated)
    class(RVAL) <- "normmix"
  }

  RVAL
}



is.normmix <- function(mix) {
  # check class of object, if TRUE, assuming it's created by normmix() and
  # subject to all its constraints.
  ifelse(inherits(mix, "normmix"), TRUE, FALSE)
}


is.intensity_surface <- function(mix) {
  # check class of object, if TRUE, assuming it's created by normmix() and
  # subject to all its constraints.
  ifelse(inherits(mix, "intensity_surface"), TRUE, FALSE)
}


#' @export
print.normmix <- function(mix) {
  cat("Normal Mixture with", mix$m, "component(s).\n")
}

#' @export
print.intensity_surface <- function(mix) {
  pre <- ifelse(mix$estimated, "Estimated", "Theoretical")
  cat(pre, "number of points over window:", mix$intensity, "\n")
  print(mix$window)
  NextMethod()
}

#' @export
print.summary_normmix <- function(mix_df) {
  with(mix_df,
       for (i in unique(comp)) {
         cat("Component", i, "is centered at", "[",
             round(mu[comp == i][1], 2), ",", round(mu[comp == i][2], 2), "]",
             "with probability", round(ps[comp == i][1], 4), "\n")
         prmatrix(round(mix_df[comp == i, 4:5], 4), rowlab = rep("", 2),
                  collab = c("covariance", "matrix:"))
       }
  )
}

#' @param mix An object of class \code{normmix}
#' @rdname normmix
#' @export
summary.normmix <- function(mix) {
  comps <- vector("list", mix$m)
  for (i in 1:mix$m) {
    comps[[i]] <- data.frame(comp = i, ps = mix$ps[[i]],
                             mu = mix$mus[[i]], sigma = mix$sigmas[[i]])
  }
  RVAL <- do.call(rbind, comps)
  class(RVAL) <- c("summary_normmix", oldClass(RVAL))
  RVAL
}

#' @rdname normmix
#' @export
summary.intensity_surface <- function(mix) {
  print(mix)
  NextMethod()
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
#' @param df Degree of freedom in generating random matrix from Wishart
#' distribution.
#' @param rand_m Whether the number of components are random or fixed (default).
#' @param xlim,ylim vector of length two, the limit are used to sample the mu's
#' from a uniform distribution.
#'
#' @return Object of class \code{normmix}.
#' @export
#' @examples
#' mix1 <- rnormmix(m = 3, sig0 = .1, df = 5)
#' summary(mix1)
#'
#' mix2 <- rnormmix(m = 5, sig0 = .1, df = 5, rand_m = TRUE, ylim = c(0, 5))
#' summary(mix2)
#'
rnormmix <- function(m, sig0, df, rand_m = FALSE,
                     xlim = c(0, 1), ylim = c(0, 1)) {
  if (rand_m) {
    # number of components is random
    m <- sample(1:m, 1)
  }

  gen_ps <- rdirichlet(1, rep(1, m))
  gen_mu <- cbind(x = runif(m, xlim[1], xlim[2]),
                  y = runif(m, ylim[1], ylim[2]))
  gen_sigma <- stats::rWishart(m, df, sig0 * diag(2))

  mus <- vector(mode = "list", length = m)
  sigmas <- vector(mode = "list", length = m)

  for (k in 1:m) {
    mus[[k]] <- gen_mu[k, ]
    sigmas[[k]] <- gen_sigma[, , k]
  }

  normmix(gen_ps, mus, sigmas)
}
#' Convert normal mixture to intensity surface
#'
#' This function can convert a normmix object into an intensity surface. It can
#' also be used to change intensity or window of an intensity_surface object.
#'
#' If the class of mix is normmix, lambda and win are used to convert mix to an
#' intensity surface class. If the class of mix is intensity_surface, lambda
#' and win are used to change the original setting of lambda and win.
#'
#' @param mix object with class normmix or intensity_surface.
#' @param lambda optional parameter of average intensity.
#' @param win optional parameter of domain.
#'
#' @examples
#' mix1 <- normmix(ps = c(.3, .7),
#'                 mus = list(c(0.2, 0.2), c(.8, .8)),
#'                 sigmas = list(.01*diag(2), .01*diag(2)))
#'
#' intsurf1 <- normmix(ps = c(.3, .7),
#'                     mus = list(c(0.2, 0.2), c(.8, .8)),
#'                     sigmas = list(.01*diag(2), .01*diag(2)),
#'                     lambda = 100,
#'                     win = square(1))
#'
#' # from normmix
#' to_int_surf(mix1, lambda = 100, win = square(1))
#'
#' # from intensity_surface
#' to_int_surf(intsurf1, win = square(2))
#' to_int_surf(intsurf1, lambda = 50)
#'
#' # to return a normmix, win will be ignored
#' to_int_surf(intsurf1, win = square(2), return_normmix = TRUE)
#' @export
to_int_surf <- function(mix, lambda = NULL, win = NULL,
                        return_normmix = FALSE) {
  intsurf <- mix

  if (is.intensity_surface(mix)) {
    # input is intensity surface
    if (!missing(lambda)) {
      stopifnot(lambda > 0)
      intsurf$intensity <- lambda
    }
    if (!missing(win)) {
      stopifnot(spatstat::is.owin(win))
      intsurf$window <- win
    }

  } else if (is.normmix(mix)) {
    if (!return_normmix) {
      # input is normmix, create intensity surface
      stopifnot(!missing(lambda) & !missing(win))
      intsurf <- normmix(mix$ps, mix$mus, mix$sigmas,
                         lambda = lambda, win = win)
    }
  } else {
    stop("mix must be of class normmix or intensity surface")
  }

  if (return_normmix) {
    intsurf <- normmix(intsurf$ps, intsurf$mus, intsurf$sigmas)
  }

  intsurf
}

