#' Generate a spatial point pattern from normal mixture.
#'
#' This function generates a point pattern from a given normal mixture or
#' intensity surface.
#'
#' If an intensity surface is passed to \code{intsurf}, it will generate mixture
#' directly. If you don't have an intensity surface beforehand, you can pass a
#' normal mixture of class \code{normmix} and specify lambda and window as
#' additianl parameters. Even if you have an intensity surface and passed it to
#' \code{rsppmix()}, you can still overwrite it's intensity and window with
#' additional parameters. See examples for details.
#'
#' The number of points \code{n} follows Poisson distribution with intensity
#' lambda * area of window.
#'
#' When \code{truncate = TRUE}, a point pattern with \code{n} points will be
#' generated from the mixture first. Then if not all the points are in the
#' domain, it will generate more points until there are exactly \code{n} points
#' in the domain. If \code{truncate = FALSE} is set, the returned point pattern
#' will not check whether the points are inside the domain.
#'
#' @param intsurf Object of class \code{intensity_surface} or \code{normmix}.
#' @param truncate Whether to truncate the points outside the domain,
#' default to TRUE.
#' @param ... Further parameters passed to \code{to_int_surf()}.
#'
#' @return A point pattern of class \code{c("sppmix", "ppp")}.
#' @export
#' @examples
#' rsppmix(intsurf = demo_intsurf)
#'
#' # overwrite lambda or win
#' rsppmix(demo_intsurf, lambda = 200)
#' rsppmix(demo_intsurf, win = square(2))
#'
#' # use normmix with additional parameters
#' rsppmix(demo_mix, lambda = 100, win = square(1))
#'
#' # turn off truncation
#' rsppmix(intsurf = demo_intsurf, truncate = FALSE)
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
      spp[[k]] <- cbind(snd, k)
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

  spp <- spp[sample(1:nrow(spp), size = nrow(spp)), ]

  RVAL <- as.ppp(spp[, 1:2], W=win, check = truncate)
  RVAL$comp <- spp[, 3]
  class(RVAL) <- c("sppmix", "ppp")
  return(RVAL)
}
