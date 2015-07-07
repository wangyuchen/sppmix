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

rnormmix <- function(m, sig0, sigdf,
                        win=spatstat::square(1),
                        rand_m=FALSE) {
  if (!spatstat::is.owin(win)) {
    stop("win must be of class owin.")
  }
  if (rand_m == TRUE) {
    # number of components is random
    m <- sample(1:m, 1)
  }

  ps <- MCMCpack::rdirichlet(1, rep(1, m))
  mus <- list()
  sigmas <- list()

  for (k in 1:m) {
    mus[[k]] <- c(runif(1, win$xrange[1], win$xrange[2]),
                runif(1, win$yrange[1], win$yrange[2]))
    sigmas[[k]] <- rWishart(1, sigdf, sig0 * diag(2))[, , 1]
  }

  return(normmix(ps, mus, sigmas))
}

dnormmix <- function(x, y, mix, win = spatstat::square(1), truncate = TRUE) {
  approx <- approx_normmix(mix, win)
  den <- matrix(NA_real_, length(x), mix$m)
  for (k in 1:mix$m) {
    # every row of den is for a point
    # every col of den is for a component
    den[, k] <- mvtnorm::dmvnorm(cbind(x, y), mix$mus[[k]], mix$sigmas[[k]])
    den[, k] <- den[, k] * mix$ps[k]
    if (truncate) {
      den[, k] <- den[, k] / approx[k]
    }
  }
  return(rowSums(den))
}

summary.normmix <- function(mix) {

}

plot.normmix <- function(mix) {

}
