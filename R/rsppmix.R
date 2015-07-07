rsppmix <- function(lambda, mix, win=spatstat::square(1), truncate=TRUE) {
  if (!is.normmix(mix)) {
    stop("mix must be an object of class normmix.")
  }

  n <- rpois(1, lambda)
  if (n == 0) {
    stop("0 points in the requested point pattern")
  }

  gen_n_from_mix <- function(n, mix) {
    comp <- sample(1:mix$m, size = n, replace=TRUE, prob=mix$ps)
    spp <- mvtnorm::rmvnorm(sum(comp == 1),
                            mix$mus[[1]], mix$sigmas[[1]])

    if (mix$m >= 1) {
      for (k in 2:mix$m) {
        samples <- mvtnorm::rmvnorm(sum(comp == k),
                                    mix$mus[[k]], mix$sigmas[[k]])
        spp <- rbind(spp, samples)
      }
    }
    return(spp)
  }

  spp <- gen_n_from_mix(n, mix)

  if (truncate == TRUE) {
    while (sum(spatstat::inside.owin(spp[, 1], spp[, 2], win)) < n) {
      spp <- rbind(spp, gen_n_from_mix(n, mix))
    }
    spp <- spp[spatstat::inside.owin(spp[, 1], spp[, 2], win), ][1:n, ]
  }

  return(as.ppp(spp, W=win, check = truncate))
}
