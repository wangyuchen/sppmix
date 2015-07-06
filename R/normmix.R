normmix <- function(ps, mus, sigmas) {
  if (length(ps) == 0) {
    stop("Mixture with 0 components")
  }

  if (sum(ps) != 1) {
    stop("Component probabilities must sum to 1.")
  }

  if (length(ps) != length(mus) | length(ps) != length(sigmas)) {
    stop("Number of components mismatch.")
  }

  RVAL <- list(num = length(ps), ps = ps, mus = mus, sigmas = sigmas)
  class(RVAL) <- "normmix"
  return(RVAL)
}

is.normmix <- function(mix) {
  # check class of object, if TRUE, assuming it's created by normmix() and
  # subject to all its constraints.
  return(ifelse(class(mix) == "normmix", TRUE, FALSE))
}

print.normmix <- function(mix) {
  print(paste("Normal Mixture with", mix$num, "component"))
}

summary.normmix <- function(mix) {

}

plot.normmix <- function(mix) {

}
