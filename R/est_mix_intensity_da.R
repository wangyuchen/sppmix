#' @export
est_mix_intensity <- function(pattern, win, m, L = 1000, burnin = 200,
                              truncate = TRUE) {
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

  sample_mu <- function(j, sigma) {
    sum1 <- count_ind(zmultinom, which = j)
    if (sum1 > 0) newmu <- colMeans(pp[zmultinom == j, ]) else newmu <- c(0, 0)

    invsig1 <- inv(sigma[, , j])
    cov1 <- inv(sum1*invsig1 + kappa)
    mu1 <- cov1 %*% (sum1*invsig1 %*% newmu + kappa %*% ksi)

    propmu <- mvtnorm::rmvnorm(1, mu1, cov1)
    return(propmu)
  }

  sample_sigma <- function(j, mu) {
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

  pp <- as.data.frame(pattern)
  n <- npoints(pattern)
  if (truncate==TRUE) {
    pattern <- pattern[spatstat::inside.owin(pattern$x,pattern$y,win)]
  }
  if (L <= burnin) {
    stop("wrong L or burnin")
  }
  MHjump <- 0
  ksi <- colMeans(pp)
  R1 <- diff(range(pp$x))
  R2 <- diff(range(pp$y))
  kappa <- diag(c(1/R1^2,1/R2^2))
  kappainv <- diag(c(R1^2,R2^2))
  a <- 3
  g <- 0.3
  gam <- 1
  hmat <- 100*g/a*kappa

  #starting value
  mus <- list(pp[sample(1:n, size = m, replace = T), ])
  sigmas <- list(array(kappainv,dim = c(2, 2, m)))
  ps <- matrix(1/m, L, m)

  #posterior sample for lambda
  alambda <- 1
  blambda <- 10000
  lambdas <- rgamma(L, n + alambda, scale = blambda/(blambda+1))
  meanlambda <- mean(lambdas)
  lognfac <- log(sqrt(2*pi*n)) + n*log(n) - n
  logdens1 <- -lognfac + n*log(meanlambda) - meanlambda

  zmultinom <- sample(1:m, size = n, replace = T)

  ## start main mcmc ##
  pb <- txtProgressBar(min = 1, max = L, initial = 2)
  el2 <- el3 <- el4 <- 0

  for (i in 2:L) {
    setTxtProgressBar(pb, i)

    #sample B matrix
    beta <- sample_beta(m, sigmas[[i-1]])

    approx <- rep(1, m)

    #sample mus and sigmas
    mus <- append(mus, list(mus[[i-1]]))
    sigmas <- append(sigmas, list(sigmas[[i-1]]))

    # sample mus
    mix_old_mu <- as.normmix(ps[i-1, ], mus[[i-1]], sigmas[[i-1]])
    approx_old_mu <- approx_normmix(mix_old_mu, win)

    propmus <- t(sapply(1:m, sample_mu, sigma = sigmas[[i-1]]))

    # truncate
    if (truncate == TRUE) {
      mix_prop_mu <- as.normmix(ps[i-1, ], propmus, sigmas[[i-1]])
      approx_prop_mu <- approx_normmix(mix_prop_mu, win)
      ratio <- approx_old_mu / approx_prop_mu
    } else {
      ratio <- 1
    }

    accept <- runif(m) < ratio
    mus[[i]][accept, ] <- propmus[accept, ]

    # sample sigmas

    propsigmas <- sapply(1:m, sample_sigma, mu = mus[[i]])

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

    accept <- runif(m) < ratio
    sigmas[[i]][, , accept] <- propsigmas[, accept]


    # sample ps
    ds <- gam + count_ind(zmultinom, total = m)
    ps[i, ] <- rdirichlet(1, ds)

    # sample zij
    mix <- as.normmix(ps[i, ], mus[[i]], sigmas[[i]])
    den <- matrix(NA_real_, n, mix$m)
    for (k in 1:mix$m) {
      den[, k] <- mvtnorm::dmvnorm(pp, mix$mus[[k]], mix$sigmas[[k]])
      den[, k] <- den[, k] * mix$ps[k]
    }

    if (truncate == TRUE) {
      consts <- approx_normmix(mix, win)
    } else {
      consts <- rep(1,m)
    }


    qij <- t(apply(den, 1, function(x) x / sum(x)))
    qij <- scale(qij,center = F, scale = consts)
    propz <- apply(qij, 1, sample, x = 1:m, size = 1, replace = T)

    ratio <- ifelse(any(count_ind(zmultinom, total = m) < 2), 0, 1)

    if (runif(1) < ratio) {
      MHjump <- MHjump + 1
      zmultinom <- propz
    } else {
      sigmas[[i]] <- sigmas[[i-1]]
      mus[[i]] <- mus[[i-1]]
      ps[i, ] <- ps[i-1, ]
    }


  }

  close(pb)

  postmus <- Reduce("+",mus[-(1:burnin)])/(L - burnin)
  postsigmas <- Reduce("+",sigmas[-(1:burnin)])/(L - burnin)
  postps <- colMeans(ps[-(1:burnin), ])
  post.mix <- as.normmix(postps, postmus, postsigmas)
  return(post.mix)

}
