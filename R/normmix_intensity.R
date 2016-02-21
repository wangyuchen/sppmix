normmix_intensity <- function(normmix, lambda) {
  if (!is.normmix(normmix))
    stop("normmix must be an object created by normmix().")

  if (lambda <= 0) stop("Intensity must be greater than 0.")

  RVAL <- c(normmix, lambda = lambda)
  class(RVAL) <- c("normmix_intensity", "normmix")
  RVAL
}
