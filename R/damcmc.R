est_mix_damcmc <- function(pattern, m, truncate = TRUE, L = 1000) {

  inv <- function(x) {
    # for faster 2 by 2 matrix inverse
    matrix(c(x[4], -x[2], -x[3], x[1]), 2, 2) / (x[1] * x[4] - x[2] * x[3])
  }

  count_ind <- function(zmultinom, total, which = 1:total) {
    # returns a vector indicating how many points are in the components
    colSums(sapply(which, function(j) zmultinom == j))
  }

  sample_beta <- function (m, sigma) {
    # internal function for sample beta
    invsigmas <- apply(sigma, 3, inv)
    sumsig <- matrix(rowSums(invsigmas), 2, 2)
    ps1 <- inv(2*hmat + 2*sumsig)
    return(stats::rWishart(1, 2*g + 2*m*a, ps1)[, , 1])
  }

  sample_mu <- function(j, sigma, zmultinom) {
    sum1 <- count_ind(zmultinom, which = j)
    if (sum1 > 0) {
      newmu <- colMeans(pp[zmultinom == j, ])
    } else {
      newmu <- c(0, 0)
    }

    invsig1 <- inv(sigma[, , j])
    cov1 <- inv(sum1*invsig1 + kappa)
    mu1 <- cov1 %*% (sum1*invsig1 %*% newmu + kappa %*% ksi)

    propmu <- mvtnorm::rmvnorm(1, mu1, cov1)
    return(propmu)
  }

  as.normmix <- function(ps, mus, sigmas) {
    # coerce mus and sigmas in DAMCMC format into normmix
    mu <- list()
    sigma <- list()
    for (i in 1:nrow(mus)) {
      mu[[i]] <- as.numeric(mus[i, ])
      sigma[[i]] <- sigmas[, , i]
    }
    return(normmix(ps, mu, sigma))
  }

  sample_sigma <- function(j, mu, beta, zmultinom) {
    # sample sigma from proposal distribution for one component
    sum1 <- sum(zmultinom == j)
    xmu <- scale(pp[zmultinom == j,], as.numeric(mu[j, ]), scale = F)

    # sumxmu is a 2 by 2 matrix
    sumxmu <- crossprod(xmu)

    ps2 <- inv(2 * beta + sumxmu)
    invsig11 <- stats::rWishart(1, 2*a + sum1, ps2)[, , 1]
    propsigma <- inv(invsig11)
    return(propsigma)
  }

  win <- pattern$window
  xlim <- win$xrange
  ylim <- win$yrange

  # Data truncation
  if (truncate==TRUE) {
    pattern <- pattern[spatstat::inside.owin(pattern$x, pattern$y, win)]
  }

  pp <- as.data.frame(pattern)
  n <- spatstat::npoints(pattern)


  ksi <- colMeans(pp)
  R1 <- diff(range(pp$x))
  R2 <- diff(range(pp$y))
  kappa <- diag(c(1/R1^2,1/R2^2))
  kappainv <- diag(c(R1^2,R2^2))
  invK <- inv(100*kappainv)
  a <- 3
  g <- 0.3
  gam <- 1
  hmat <- 100*g/a*kappa
  loglik <- rep(0, L)


  #starting value
  betas <- sigmas <- mus <- vector("list", length = L)

  mus[[1]] <- pp[sample(1:n, size = m, replace = T), ]
  sigmas[[1]] <- array(kappainv, dim = c(2, 2, m))
  betas[[1]] <- matrix(0, 2, 2)
  ps <- matrix(1/m, L, m)

  #posterior sample for lambda
  alambda <- 1
  blambda <- 10000
  lambdas <- rgamma(L, n + alambda, scale = blambda/(blambda+1))

  zmultinom <- sample(1:m, size = n, replace = T)
  propz <- zmultinom
  qij <- matrix(1/m, nrow = n, ncol = m)

  ## start main mcmc ##
  pb <- txtProgressBar(min = 1, max = L, initial = 2, style = 3)

  time1 <- time2 <- time3 <- time4 <- 0

  for (i in 2:L) {
    t1 <- Sys.time()

    setTxtProgressBar(pb, i)

    #sample B matrix
    betas[[i]] <- sample_beta(m, sigmas[[i - 1]])

    approx <- rep(1, m)

    #sample mus and sigmas
    mus[[i]] <- mus[[i - 1]]
    sigmas[[i]] <- sigmas[[i - 1]]


    # sample mus
    propmus <- t(sapply(1:m, sample_mu, sigma = sigmas[[i - 1]],
                        zmultinom = zmultinom))


    # truncate
    if (truncate == TRUE) {
      mix_old_mu <- as.normmix(ps[i - 1, ], mus[[i - 1]], sigmas[[i - 1]])
      approx_old_mu <- approx_normmix(mix_old_mu, xlim, ylim)

      mix_prop_mu <- as.normmix(ps[i - 1, ], propmus, sigmas[[i - 1]])
      approx_prop_mu <- approx_normmix(mix_prop_mu, xlim, ylim)

      ratio <- approx_old_mu / approx_prop_mu
    } else {
      ratio <- 1
    }

    accept <- runif(m) < ratio
    if (any(accept == TRUE)) {
      mus[[i]][accept, ] <- propmus[accept, ]
    }

    t2 <- Sys.time()
    time1 <- time1 + t2 - t1

    # sample sigmas
    propsigmas <- sapply(1:m, sample_sigma,
                         mu = mus[[i]], betas[[i]], zmultinom = zmultinom)

    if (truncate == TRUE) {
      if (all(accept)) {
        # all mus proposed are accepted
        approx_old_sigma <- approx_prop_mu
      } else {
        mix_old_sigma <- as.normmix(ps[i-1, ], mus[[i]], sigmas[[i-1]])
        approx_old_sigma <- approx_normmix(mix_old_sigma, xlim, ylim)
      }

      mix_prop_sigma <- as.normmix(ps[i-1, ], mus[[i]],
                                   array(propsigmas, dim = c(2, 2, m)))
      approx_prop_sigma <- approx_normmix(mix_prop_sigma, xlim, ylim)

      ratio <- approx_old_sigma / approx_prop_sigma
    } else {
      ratio <- 1
    }

    accept <- runif(m) < ratio
    sigmas[[i]][, , accept] <- propsigmas[, accept]

    t3 <- Sys.time()
    time2 <- time2 + t3 - t2

    # sample ps
    ds <- gam + count_ind(zmultinom, total = m)
    ps[i, ] <- rdirichlet(1, ds)

    # sample zij
    mix <- as.normmix(ps[i, ], mus[[i]], sigmas[[i]])
    approx <- approx_normmix(mix, xlim, ylim)
    den <- matrix(NA_real_, n, m)
    for (k in 1:mix$m) {
      den[, k] <- mvtnorm::dmvnorm(pp, mix$mus[[k]], mix$sigmas[[k]])
      den[, k] <- den[, k] * mix$ps[k]
      if (truncate) {
        den[, k] <- den[, k] / approx[k]
      }
    }

    t4 <- Sys.time()
    time3 <- time3 + t4 - t3

    # qij is n by m
    if (m == 1){
      qij <- matrix(apply(den, 1, function(x) x / sum(x)), n, 1)
    } else{
      qij <- t(apply(den, 1, function(x) x / sum(x)))
    }

    propz <- apply(qij, 1, sample, x = 1:m, size = 1, replace = T)

    if (all(count_ind(propz, total = m) >= 2)) {
      # accept
      zmultinom <- propz
    } else {
      # reject
      ps[i, ] <- ps[i-1, ]
      sigmas[[i]] <- sigmas[[i-1]]
      mus[[i]] <- mus[[i-1]]
    }


    t5 <- Sys.time()
    time4 <- time4 + t5 - t4

  }

  close(pb)

  print(c(time1, time2, time3, time4))

  RVAL <- list(lambdas = lambdas,
               ps = ps,
               mus = mus,
               sigmas = sigmas)
  return(invisible(RVAL))
}
