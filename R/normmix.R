#' Two-dimensional mixture with normal components.
#'
#' Constructor function of \code{normmix} class. Create a mixture in
#' two-dimensional with bivariate normal components. If additional parameters
#' lambda and window are set, it will create an intensity surface.
#'
#' @param ps Vector of component probabilities.
#' @param mus A list where every element is a vector of length 2, defining the
#' center of each component.
#' @param sigmas A list where every element is a 2 by 2 covariance matrix,
#' defining the covariance for each component.
#' @param lambda Optional parameter of average number of points. If set with
#' window, the returned object will be an intensity surface.
#' @param win Optional parameter of window. Must be set with intensity value to
#' create an intensity surface.
#' @param estimated Whether it's an estimated mixture. By default it's set to
#' FALSE. When using an estimated mixture, it should be set to TRUE.
#' @param ... Currently omitted.
#'
#' @return An object of class "normmix" containing the following components:
#'  \item{m}{Number of components.}
#'  \item{ps}{Vector of component probabilities.}
#'  \item{mus}{List of mean vectors of components.}
#'  \item{sigmas}{List of covariance matrix of components.}
#'  \item{intensity}{Optional intensity value if lambda is provided when calling
#'  \code{normmix}.}
#'  \item{window}{Optional window object of class \code{\link[spatstat]{owin}}.}
#'  \item{estimated}{Whether the normal mixture is estimated.}
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
  inherits(mix, "normmix")
}


is.intensity_surface <- function(mix) {
  inherits(mix, "intensity_surface")
}


#' @export
print.normmix <- function(x, ...) {
  cat("Normal Mixture with", x$m, "component(s).\n")
}

#' @export
print.intensity_surface <- function(x, ...) {
  pre <- ifelse(x$estimated, "Estimated average", "Expected")
  cat(pre, "number of points over window:", x$intensity, "\n")
  print(x$window)
  NextMethod()
}

#' @export
print.summary_normmix <- function(x, ...) {
  with(x,
       for (i in unique(comp)) {
         cat("Component", i, "is centered at", "[",
             round(mu[comp == i][1], 2), ",", round(mu[comp == i][2], 2), "]",
             "with probability", round(ps[comp == i][1], 4), "\n")
         prmatrix(round(x[comp == i, 4:5], 4), rowlab = rep("", 2),
                  collab = c("covariance", "matrix:"))
       }
  )
}

#' @param object An object of class \code{normmix}
#' @rdname normmix
#' @export
summary.normmix <- function(object, ...) {
  comps <- vector("list", object$m)
  for (i in 1:object$m) {
    comps[[i]] <- data.frame(comp = i, ps = object$ps[[i]],
                             mu = object$mus[[i]], sigma = object$sigmas[[i]])
  }
  RVAL <- do.call(rbind, comps)
  class(RVAL) <- c("summary_normmix", oldClass(RVAL))
  RVAL
}

#' @rdname normmix
#' @export
summary.intensity_surface <- function(object, ...) {
  print(object)
  NextMethod()
}


#' Convert normal mixture to intensity surface.
#'
#' This function can convert a normmix object into an intensity surface. It can
#' also be used to change intensity or window of an intensity_surface object.
#'
#' If the class of mix is normmix, lambda and win are used to convert mix to an
#' intensity surface class. If the class of mix is intensity_surface, lambda
#' and win are used to change the original setting of lambda and win.
#'
#' @param mix Object with class normmix or intensity_surface.
#' @param lambda Optional parameter of average intensity.
#' @param win Optional parameter of domain.
#' @param return_normmix Whether to return a normal mixture. (disgard lambda
#' and win).
#'
#' @examples
#' # from normmix
#' to_int_surf(demo_mix, lambda = 100, win = square(1))
#'
#' # from intensity_surface
#' to_int_surf(demo_intsurf, win = square(2))
#' to_int_surf(demo_intsurf, lambda = 50)
#'
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
                         lambda = lambda, win = win, estimated = mix$estimated)
    }
  } else {
    stop("mix must be of class normmix or intensity surface")
  }

  if (return_normmix) {
    intsurf <- normmix(intsurf$ps, intsurf$mus, intsurf$sigmas,
                       estimated = intsurf$estimated)
  }

  intsurf
}


#' Demo objects
#'
#' Demo objects from the classes provided by this package.
#'
#' @examples
#' demo_mix <- normmix(ps = c(.3, .7), mus = list(c(0.2, 0.2), c(.8, .8)),
#'                     sigmas = list(.01*diag(2), .01*diag(2)))
"demo_mix"

#' @rdname demo_mix
#' @examples
#' demo_intsurf <- normmix(ps = c(.3, .7), mus = list(c(0.2, 0.2), c(.8, .8)),
#'                         sigmas = list(.01*diag(2), .01*diag(2)),
#'                         lambda = 100, win = square(1))
"demo_intsurf"





