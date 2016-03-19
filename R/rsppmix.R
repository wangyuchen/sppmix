#' Generate a spatial point pattern from normal mixture.
#'
#' This function generates a point pattern from a given intensity surface.
#'
#'
#' This function generates a spatial point pattern from a intensity surface.
#' Optionally, it can generate point pattern from a normal mixture if lambda
#' and window are specified.
#' The number of points \code{n} follows Poisson distribution with intensity
#' lambda * area of window.
#'
#' When \code{truncate = TRUE}, a point pattern with \code{n} points will be
#' generated from the mixture first. Then if not all the points are in the
#' domain, it will generate another \code{n} points until there are more than
#' \code{n} points in the domain. The first \code{n} points are returned as
#' the generated spatial point pattern.
#'
#' If \code{truncate = FALSE} is set, the returned point pattern will not check
#' whether the points are inside the domain.
#'
#' @param intsurf Object of class intensity_surface
#' @param truncate Whether to truncate the points outside the domain,
#' default to TRUE.
#' @param ... Further parameters passed to \code{to_int_surf()}.
#'
#' @return A point pattern of class \code{c("sppmix", "ppp")}.
#' @export
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
#' rsppmix(intsurf = intsurf1)
#'
#' # overwrite lambda or win
#' rsppmix(intsurf1, lambda = 200)
#' rsppmix(intsurf1, win = square(2))
#'
#' # use normmix with additional parameters
#' rsppmix(mix1, lambda = 100, win = square(1))
#'
#' # turn off truncation
#' rsppmix(intsurf = intsurf1, truncation = FALSE)
#'
rsppmix <- function(intsurf, truncate = TRUE, ...) {

  intsurf <- to_int_surf(intsurf, ...)
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
