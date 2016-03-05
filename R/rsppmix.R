#' Generate a spatial point pattern from normal mixture.
#'
#' Generate a spatial point pattern from a given mixture with normal components.
#'
#' @param lambda Intensity of the spatial point pattern
#' @param intsurf Object of class intensity_surface
#' @param win Object of class spatstat::owin
#' @param truncate Whether to truncate the points outside the domain,
#' default to TRUE.
#' @details
#' If \code{truncate = FALSE} is set, the returned point pattern will not check
#' whether the points are inside the domain.
#'
#' The number of points \code{n} follows Poisson distribution with intensity
#' lambda * area of window.
#'
#' When \code{truncate = TRUE}, a point pattern with \code{n} points will be
#' generated from the mixture first. Then if not all the points are in the
#' domain, it will generate another \code{n} points until there are more than
#'  \code{n} points in the domain. The first \code{n} points are returned as
#'   the generated spatial point pattern.
#'
#' @return A point pattern of class \code{c("sppmix", "ppp")}.
#' @export
#' @examples
#' # generate a intensity surface
#' intsurf1 <- normmix(ps=c(.5, .5),
#'                 mus=list(c(.2, .2), c(.7, .7)),
#'                 sigmas=list(.01*diag(2), .02*diag(2)), 200, square(1))
#' rsppmix(intsurf1)
#'
#' # If truncate = FALSE it will generate points outside the window
#' rsppmix(intsurf1, truncate = FALSE)


rsppmix <- function(intsurf, truncate = TRUE, lambda, win) {

  if (!is.intensity_surface(intsurf)) {
    if (!is.normmix(intsurf)) stop("mix must be an object of class normmix.")
    stopifnot(missing(lambda) | missing(win))
    intsurf <- to_int_surf(intsurf, lambda, win)
  }

  win <- intsurf$window

  n <- rpois(1, intsurf$intensity)

  if (n == 0) {
    stop("Intensity value too small.")
  }

  gen_n_from_mix <- function(n, mix) {
    comp <- sample(1:mix$m, size = n, replace = TRUE, prob = mix$ps)
    spp <- vector("list", length(unique(comp)))
    for (k in unique(comp)) {
      snd <- mvtnorm::rmvnorm(sum(comp == k),
                                   mix$mus[[k]], mix$sigmas[[k]])
      spp[[k]] <- cbind(snd, comp = k)
    }
    return(do.call(rbind, spp))
  }

  spp <- gen_n_from_mix(n, intsurf)

  if (truncate == TRUE) {
    while (sum(spatstat::inside.owin(spp[, 1], spp[, 2], win)) < n) {
      spp <- rbind(spp, gen_n_from_mix(n, intsurf))
    }
    spp <- spp[spatstat::inside.owin(spp[, 1], spp[, 2], win), ][1:n, ]
  } else {
    warning(paste(sum(!spatstat::inside.owin(spp[, 1], spp[, 2], win)),
                  "points are outside window."))
  }

  RVAL <- as.ppp(spp[, 1:2], W=win, check = truncate)
  RVAL$comp <- spp[, 3]
  class(RVAL) <- c("sppmix", "ppp")
  return(RVAL)
}
