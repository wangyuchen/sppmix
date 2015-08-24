#' @export
est_mix_intensity <- function(pattern, win, m, L = 1000, burnin = 200,
                              truncate = TRUE, marginal = FALSE) {
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
    return(rWishart(1, 2*g + 2*m*a, ps1)[, , 1])
  }

  den_beta <- function (m, beta, sigma) {
    # internal function for sample beta
    invsigmas <- apply(sigma, 3, inv)
    sumsig <- matrix(rowSums(invsigmas), 2, 2)
    ps1 <- inv(2*hmat + 2*sumsig)
    return(MCMCpack::dwish(beta, 2*g + 2*m*a, ps1))
  }

  sample_mu <- function(j, sigma) {
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

  sample_sigma <- function(j, mu, beta = beta) {
    # sample sigma from proposal distribution for one component
    sum1 <- sum(zmultinom == j)
    xmu <- scale(pp[zmultinom == j,], as.numeric(mu[j, ]), scale = F)

    # sumxmu is a 2 by 2 matrix
    sumxmu <- crossprod(xmu)

    ps2 <- inv(2 * beta + sumxmu)
    invsig11 <- rWishart(1, 2*a + sum1, ps2)[, , 1]
    propsigma <- inv(invsig11)
    return(propsigma)
  }

    # estimated the posterior density of sigmas
  den_sigma <- function(j, mu, beta, sigma) {
    # sample sigma from proposal distribution for one component
    nj <- sum(zmultinom == j)
    xmu <- scale(pp[zmultinom == j,], as.numeric(mu[j, ]), scale = F)

    # sumxmu is a 2 by 2 matrix
    sumxmu <- crossprod(xmu)

    ps2 <- inv(2 * beta + sumxmu)
   den <- MCMCpack::dwish(inv(sigma[, , j]), 2*a + nj, ps2)
    return(den)
  }

  # Data truncation
  if (truncate==TRUE) {
    pattern <- pattern[spatstat::inside.owin(pattern$x, pattern$y, win)]
  }

  pp <- as.data.frame(pattern)
  n <- npoints(pattern)

  if (L <= burnin) {
    stop("wrong L or burnin")
  }

  MHjump <- 0
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
  perms <-  matrix(unlist(permn(m)), ncol = m, byrow = T)
  nperms <-  dim(perms)[1]
  loglik <- rep(0, L)

  #starting value
  mus <- list(pp[sample(1:n, size = m, replace = T), ])
  sigmas <- list(array(kappainv, dim = c(2, 2, m)))
  betas <- list(matrix(0,2,2))
  ps <- matrix(1/m, L, m)
  denosum <- 0

  #posterior sample for lambda
  alambda <- 1
  blambda <- 10000
  lambdas <- rgamma(L, n + alambda, scale = blambda/(blambda+1))
  meanlambda <- mean(lambdas)

  lognfac <- log(sqrt(2*pi*n)) + n*log(n) - n
  logdens1 <- -lognfac + n*log(meanlambda) - meanlambda

  zmultinom <- sample(1:m, size = n, replace = T)
  propz <- zmultinom
  qij <- matrix(1/m, nrow = n, ncol = m)

  ## start main mcmc ##
  pb <- txtProgressBar(min = 1, max = L, initial = 2, style = 3)
  el2 <- el3 <- el4 <- 0

  for (i in 2:L) {
    setTxtProgressBar(pb, i)

    #sample B matrix
    beta <- sample_beta(m, sigmas[[i-1]])
    betas <- append(betas, list(beta))

    approx <- rep(1, m)

    #sample mus and sigmas
    mus <- append(mus, list(mus[[i-1]]))
    sigmas <- append(sigmas, list(sigmas[[i-1]]))

    # sample mus
    propmus <- t(sapply(1:m, sample_mu, sigma = sigmas[[i-1]]))

    # truncate
    if (truncate == TRUE) {
      mix_old_mu <- as.normmix(ps[i-1, ], mus[[i-1]], sigmas[[i-1]])
      approx_old_mu <- approx_normmix(mix_old_mu, win)

      mix_prop_mu <- as.normmix(ps[i-1, ], propmus, sigmas[[i-1]])
      approx_prop_mu <- approx_normmix(mix_prop_mu, win)

      ratio <- approx_old_mu / approx_prop_mu
    } else {
      ratio <- 1
    }
    accept <- rep(TRUE, m)
    accept[is.nan(ratio)] <- FALSE
    accept[!is.nan(ratio)] <- (runif(m) < ratio)[!is.nan(ratio)]
    if (any(accept != FALSE)){
      mus[[i]][accept, ] <- propmus[accept, ]
    }

    # sample sigmas
    propsigmas <- sapply(1:m, sample_sigma, mu = mus[[i]], beta)

    if (truncate == TRUE) {
      if (all(accept)) {
        # all mus proposed are accepted
        approx_old_sigma <- approx_prop_mu
      } else {
        mix_old_sigma <- as.normmix(ps[i-1, ], mus[[i]], sigmas[[i-1]])
        approx_old_sigma <- approx_normmix(mix_old_sigma, win)
      }

      mix_prop_sigma <- as.normmix(ps[i-1, ], mus[[i]],
                                   array(propsigmas, dim = c(2, 2, m)))
      approx_prop_sigma <- approx_normmix(mix_prop_sigma, win)

      ratio <- approx_old_sigma / approx_prop_sigma
    } else {
      ratio <- 1
    }

    accept[is.nan(ratio)] <- FALSE
    accept[!is.nan(ratio)] <- (runif(m) < ratio)[!is.nan(ratio)]
    sigmas[[i]][, , accept] <- propsigmas[, accept]

    # sample ps
    ds <- gam + count_ind(zmultinom, total = m)
    ps[i, ] <- rdirichlet(1, ds)

    # sample zij
    mix <- as.normmix(ps[i, ], mus[[i]], sigmas[[i]])
    approx <- approx_normmix(mix, win)
    den <- matrix(NA_real_, n, m)
    for (k in 1:mix$m) {
      den[, k] <- mvtnorm::dmvnorm(pp, mix$mus[[k]], mix$sigmas[[k]])
      den[, k] <- den[, k] * mix$ps[k]
      if (truncate) {
        den[, k] <- den[, k] / approx[k]
      }
    }

    # qij is n by m
    if (m == 1){
      qij <- matrix(apply(den, 1, function(x) x / sum(x)), n, 1)
    } else{
      qij <- t(apply(den, 1, function(x) x / sum(x)))
    }
    cond <- abs(rowSums(qij) - 1) < .0001
    propz[cond] <- apply(as.matrix(qij[cond, ]), 1, sample,
                         x = 1:m, size = 1, replace = T)


    ratio <- ifelse(any(count_ind(propz, total = m) < 2), 0, 1)

    if (runif(1) < ratio) {
      if (i > burnin) {
        MHjump <- MHjump + 1
      }
      zmultinom <- propz
    } else {
      ps[i, ] <- ps[i-1, ]
      sigmas[[i]] <- sigmas[[i-1]]
      mus[[i]] <- mus[[i-1]]
    }

  # calculate loglikehood if necessary
    if (marginal == TRUE){
      loglik[i] <- sum(log(apply(den, 1, sum)))
    }
  }
  close(pb)

  postmus <- Reduce("+", mus[-(1:burnin)])/(L - burnin)
  postsigmas <- Reduce("+", sigmas[-(1:burnin)])/(L - burnin)
  postps <- colMeans(as.matrix(ps[-(1:burnin), ]))
  post_mix <- as.normmix(postps, postmus, postsigmas)

  #calculate the marginal distribution
  if (marginal == TRUE) {
    # find MLE of each parameters
  logmle_index <- order(loglik[-(1:burnin)])[(L-burnin)]
  mlemus <- mus[[logmle_index]]
  mlesigs <- sigmas[[logmle_index]]
  invK <- inv(100*kappainv)
  Bhat <- betas[[logmle_index]]
  mleps = ps[-(1:burnin), ][logmle_index, ]
  logmle1 <- loglik[-(1:burnin)][logmle_index]
  logmle2 <- sum(mvtnorm::dmvnorm(mlemus, mean = c(mean(pp[, 1]), mean(pp[, 2])),
                     sigma =invK ,log = TRUE)) +
                   sum(log(apply(mlesigs, 3, function(W) MCMCpack::dwish(W, v = 2*a,
                                                                     S = Bhat))))
                   + log(MCMCpack::dwish(W = Bhat, v = 3, S = inv(2*hmat))) +
                   +lgamma(n)
  # Posterior dentsity of Beta
  density_beta <- mean(sapply(1:(L-burnin),
                          function(x) den_beta(m = m,
                                                   beta = Bhat,
                                                   sigma = sigmas[-(1:burnin)][[x]])))
  # Run a reduced Gibbs for sigmas by given Bhat

  mus2 <- list(pp[sample(1:n, size = m, replace = T), ])
  sigmas2 <- list(array(kappainv, dim = c(2, 2, m)))
  ps2 <- matrix(1/m, L, m)
  zmultinom2 <- sample(1:m, size = n, replace = T)
  propz <- zmultinom2

  for (i in 2:L) {
    approx <- rep(1, m)

    #sample mus and sigmas
    mus2 <- append(mus2, list(mus2[[i-1]]))
    sigmas2 <- append(sigmas2, list(sigmas2[[i-1]]))

    # sample mus
    propmus <- t(sapply(1:m, sample_mu, sigma = sigmas2[[i-1]]))

    # truncate
    if (truncate == TRUE) {
      mix_old_mu <- as.normmix(ps2[i-1, ], mus2[[i-1]], sigmas2[[i-1]])
      approx_old_mu <- approx_normmix(mix_old_mu, win)

      mix_prop_mu <- as.normmix(ps2[i-1, ], propmus, sigmas2[[i-1]])
      approx_prop_mu <- approx_normmix(mix_prop_mu, win)

      ratio <- approx_old_mu / approx_prop_mu
    } else {
      ratio <- 1
    }
    accept <- rep(TRUE, m)
    accept[is.nan(ratio)] <- FALSE
    accept[!is.nan(ratio)] <- (runif(m) < ratio)[!is.nan(ratio)]
    if (any(accept != FALSE)){
      mus2[[i]][accept, ] <- propmus[accept, ]
    }

    # sample sigmas
    propsigmas <- sapply(1:m, sample_sigma, mu = mus2[[i]], beta = Bhat)

    if (truncate == TRUE) {
      if (all(accept)) {
        # all mus proposed are accepted
        approx_old_sigma <- approx_prop_mu
      } else {
        mix_old_sigma <- as.normmix(ps2[i-1, ], mus2[[i]], sigmas2[[i-1]])
        approx_old_sigma <- approx_normmix(mix_old_sigma, win)
      }

      mix_prop_sigma <- as.normmix(ps2[i-1, ], mus2[[i]],
                                   array(propsigmas, dim = c(2, 2, m)))
      approx_prop_sigma <- approx_normmix(mix_prop_sigma, win)

      ratio <- approx_old_sigma / approx_prop_sigma
    } else {
      ratio <- 1
    }

    accept[is.nan(ratio)] <- FALSE
    accept[!is.nan(ratio)] <- (runif(m) < ratio)[!is.nan(ratio)]
    sigmas2[[i]][, , accept] <- propsigmas[, accept]

    # sample ps
    ds <- gam + count_ind(zmultinom2, total = m)
    ps2[i, ] <- rdirichlet(1, ds)

    # sample zij
    mix <- as.normmix(ps2[i, ], mus2[[i]], sigmas2[[i]])
    approx <- approx_normmix(mix, win)
    den <- matrix(NA_real_, n, m)
    for (k in 1:mix$m) {
      den[, k] <- mvtnorm::dmvnorm(pp, mix$mus[[k]], mix$sigmas[[k]])
      den[, k] <- den[, k] * mix$ps[k]
      if (truncate) {
        den[, k] <- den[, k] / approx[k]
      }
    }

    # qij is n by m
    if (m == 1){
      qij <- matrix(apply(den, 1, function(x) x / sum(x)), n, 1)
    } else{
      qij <- t(apply(den, 1, function(x) x / sum(x)))
    }
    cond <- abs(rowSums(qij) - 1) < .0001
    propz[cond] <- apply(as.matrix(qij[cond, ]), 1, sample,
                         x = 1:m, size = 1, replace = T)


    ratio <- ifelse(any(count_ind(propz, total = m) < 2), 0, 1)

    if (runif(1) < ratio) {
      if (i > burnin) {
        MHjump <- MHjump + 1
      }
      zmultinom2 <- propz
    } else {
      ps2[i, ] <- ps2[i-1, ]
      sigmas2[[i]] <- sigmas2[[i-1]]
      mus2[[i]] <- mus2[[i-1]]
    }
  }
  sigm <- function(m, mu, beta, sigma) {
     logsum<- exp(sum(log(sapply(1:m, den_sigma, mu = mu,
           beta = beta, sigma = sigma))))
    return(logsum)
  }
  density_sigma <- mean(sapply(1:(L-burnin),
                               function(x) sigm(m = m,
                                                mu = mus2[-(1:burnin)][[x]],
                                                beta = Bhat,
                                                sigma = mlesigs)))

  # Run a subsequence of reduced Gibbs for mlemus by given Bhat and mlesigs

  mus3 <- list(pp[sample(1:n, size = m, replace = T), ])
  ps3 <- matrix(1/m, L, m)
  zmultinom3 <- matrix(0, L, n)
  zmultinom3[1, ] <- sample(1:m, size = n, replace = T)
  propz <- zmultinom3[1, ]

  for (i in 2:L) {
    approx <- rep(1, m)

    #sample mus
    mus3 <- append(mus3, list(mus3[[i-1]]))
    propmus <- t(sapply(1:m, sample_mu, sigma = mlesigs))

    # truncate
    if (truncate == TRUE) {
      mix_old_mu <- as.normmix(ps3[i-1, ], mus3[[i-1]], mlesigs)
      approx_old_mu <- approx_normmix(mix_old_mu, win)

      mix_prop_mu <- as.normmix(ps3[i-1, ], propmus, mlesigs)
      approx_prop_mu <- approx_normmix(mix_prop_mu, win)

      ratio <- approx_old_mu / approx_prop_mu
    } else {
      ratio <- 1
    }
    accept <- rep(TRUE, m)
    accept[is.nan(ratio)] <- FALSE
    accept[!is.nan(ratio)] <- (runif(m) < ratio)[!is.nan(ratio)]
    if (any(accept != FALSE)){
      mus3[[i]][accept, ] <- propmus[accept, ]
    }

    # sample ps
    ds <- gam + count_ind(zmultinom3[i-1, ], total = m)
    ps3[i, ] <- rdirichlet(1, ds)

    # sample zij
    mix <- as.normmix(ps3[i, ], mus3[[i]], mlesigs)
    approx <- approx_normmix(mix, win)
    den <- matrix(NA_real_, n, m)
    for (k in 1:mix$m) {
      den[, k] <- mvtnorm::dmvnorm(pp, mix$mus[[k]], mix$sigmas[[k]])
      den[, k] <- den[, k] * mix$ps[k]
      if (truncate) {
        den[, k] <- den[, k] / approx[k]
      }
    }

    # qij is n by m
    if (m == 1){
      qij <- matrix(apply(den, 1, function(x) x / sum(x)), n, 1)
    } else{
      qij <- t(apply(den, 1, function(x) x / sum(x)))
    }
    cond <- abs(rowSums(qij) - 1) < .0001
    propz[cond] <- apply(as.matrix(qij[cond, ]), 1, sample,
                         x = 1:m, size = 1, replace = T)


    ratio <- ifelse(any(count_ind(propz, total = m) < 2), 0, 1)

    if (runif(1) < ratio) {
      if (i > burnin) {
        MHjump <- MHjump + 1
      }
      zmultinom3[i, ] <- propz
    } else {
      ps3[i, ] <- ps3[i-1, ]
      mus3[[i]] <- mus3[[i-1]]
    }
  }
  mum <- function(m, mu, sigma, zmultinom) {
    invsigmas <- apply(sigma, 3, inv)
    sum1 <- count_ind(zmultinom, total = m)
    if (all(sum1 > 0)) {
      newmu <- aggregate(data.frame(x=pp[,1],y=pp[,2]), list(zmultinom), mean)
      newmu <- as.matrix(cbind(newmu$x,newmu$y))
    } else {
      newmu <- matrix(0, m, 2)
    }

    par1 <- sapply(1:m,function(x) sum1[x]*invsigmas[, x])
    par2 <- array(apply(par1, 2, function(x) matrix(x,2,2) + kappa),
                  dim = c(2, 2, m))
    cov1 <- apply(par2,3,inv)
    cov1 <- array(apply(cov1,2,function(x) matrix(x, 2, 2)), dim = c(2, 2, m))
    mu1 <- t(sapply(1:m, function(x) cov1[, , x] %*%
                      (sum1[x]*matrix(invsigmas[,x], 2, 2) %*%
                         newmu[x,] + kappa %*% ksi)))
    rval <- exp(sum(sapply(1:m, function(x) mvtnorm::dmvnorm(mu[x, ],
                       mean = mu1[x, ], sigma = cov1[, , x], log = TRUE))))
    return(rval)
  }


  density_mu <- mean(sapply(1:(L-burnin),
                            function(x) mum(m = m,
                                            mu = mlemus,
                                            sigma = mlesigs,
                                            zmultinom =
                                              zmultinom3[-(1:burnin), ][x, ])))
  # Run a sub-sequence of reduced Gibbs for mleps by given other parameters

  ps4 <- matrix(1/m, L, m)
  zmultinom4 <- matrix(0, L, n)
  zmultinom4[1, ] <- sample(1:m, size = n, replace = T)
  propz <- zmultinom4[1, ]

  for (i in 2:L) {
    approx <- rep(1, m)

    # sample ps
    ds <- gam + count_ind(zmultinom4[i-1, ], total = m)
    ps4[i, ] <- rdirichlet(1, ds)

    # sample zij
    mix <- as.normmix(ps4[i, ], mlemus, mlesigs)
    approx <- approx_normmix(mix, win)
    den <- matrix(NA_real_, n, m)
    for (k in 1:mix$m) {
      den[, k] <- mvtnorm::dmvnorm(pp, mix$mus[[k]], mix$sigmas[[k]])
      den[, k] <- den[, k] * mix$ps[k]
      if (truncate) {
        den[, k] <- den[, k] / approx[k]
      }
    }

    # qij is n by m
    if (m == 1){
      qij <- matrix(apply(den, 1, function(x) x / sum(x)), n, 1)
    } else{
      qij <- t(apply(den, 1, function(x) x / sum(x)))
    }
    cond <- abs(rowSums(qij) - 1) < .0001
    propz[cond] <- apply(as.matrix(qij[cond, ]), 1, sample,
                         x = 1:m, size = 1, replace = T)


    ratio <- ifelse(any(count_ind(propz, total = m) < 2), 0, 1)

    if (runif(1) < ratio) {
      if (i > burnin) {
        MHjump <- MHjump + 1
      }
      zmultinom4[i, ] <- propz
    } else {
      ps4[i, ] <- ps4[i-1, ]
    }
  }
  den_ps <- function(ps, zmultinom) {
    sum1 <- count_ind(zmultinom, total = m)
    den <- exp(sum((sum1 -1)*log(ps)) +
      lgamma(n + m) - sum(lgamma(sum1 + 1)))
    return(den)
  }
  density_ps <- mean(sapply(1:(L-burnin),
                       function(x) den_ps(ps4[-(1:burnin), ][x, ],
                                          zmultinom4[-(1:burnin), ][x, ])))
  # estimate the posterior density ordinate
  marginal <- exp(logmle1 + logmle2)/(density_ps*density_mu*density_sigma*
                                      density_beta)
  }
  RVAL <- list(lambda = meanlambda,
               ps = ps[-(1:burnin), ],
               mus = mus[-(1:burnin)],
               sigmas = sigmas[-(1:burnin)],
               post_mix = post_mix,
               marginal = marginal,
               accept_rate = MHjump / (L - burnin))
  class(RVAL) <- "dares"
  return(RVAL)
}
