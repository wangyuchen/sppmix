#' @export
est_mix_intensity <- function(pattern, win, m, L = 10000, burnin = 2000,
                              truncate = TRUE) {

  sample_beta <- function (m, sigma, hmat, g, a, n) {
    invsigmas <- array(dim = c(2, 2, m))
    for (j in 1:m) {
      invsigmas[, , j] <- solve(sigma[, , j])
    }
    sumsig <- apply(invsigmas, 1:2, sum)
    ps1 <- solve(2*hmat + 2*sumsig)
    return(rWishart(1, 2*g + 2*m*a, ps1))
  }

  sample_mu <- function(j, old_mu, approx_old_mu, sigma) {
    # sample mu from proposal for one component
    sum1 <- sum(zmultinom == j)
    if (sum1 > 0) {
      newmu <- colMeans(pp[zmultinom == j, ])
    } else {
      newmu <- c(0, 0)
    }

    #sample mus
    invsig1 <- solve(sigma)
    cov1 <- solve(sum1*invsig1 + kappa)
    mu1 <- cov1 %*% (sum1*invsig1 %*% newmu + kappa %*% ksi)
    propmu <- mvtnorm::rmvnorm(1, mu1, cov1)
    if (truncate == TRUE) {
      mix_prop_mu <- as.normmix(ps[i-1, ], propmu, sigmas[[i-1]])
      approx_prop_mu <- approx_normmix(mix_prop_mu, win)[j]
      ratio <- approx_old_mu / approx_prop_mu
    } else {
      ratio <- 1
    }

    if (runif(1) < ratio) new_mu <- propmu else new_mu <- old_mu
    return(new_mu)
  }

  sample_sigma <- function(j, mu, old_sigma) {
    # sample sigma from proposal distribution for one component
    sum1 <- sum(zmultinom == j)
    xmu <- as.matrix(pp[zmultinom == j,] - mu)

    # sumxmu is a 2 by 2 matrix
    sumxmu <- crossprod(as.matrix(xmu))

    ps2 <- solve(2*beta[ , , 1] + sumxmu)
    invsig11 <- rWishart(1, 2*a + sum1, ps2)
    propsigma <- solve(invsig11[, , 1])

    if (truncate == TRUE) {
      mix_old_sigma <- as.normmix(ps[i-1, ], mu, old_sigma)
      approx_old_sigma <- approx_normmix(mix_old_sigma, win)[j]
      mix_prop_sigma <- as.normmix(ps[i-1, ], mu, propsigma)
      approx_prop_sigma <- approx_normmix(mix_prop_sigma, win)[j]
      ratio <- mix_old_sigma / mix_prop_sigma
    } else {
      ratio <- 1
    }

    if (runif(1) < ratio) new_sigma <- propsigma else new_sigma <- old_sigma
    return(new_sigma)
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
  consts <- matrix(0, L, m)

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

  for (i in 2:L) {
    setTxtProgressBar(pb, i)
    #sample B matrix

    beta <- sample_beta(m, sigmas[[i-1]], hmat, g, a, n)

    approx <- rep(1, m)

    #sample mus and sigmas
    mus <- append(mus, list(matrix(NA, m, 2)))
    sigmas <- append(sigmas, list(array(NA, dim = c(2, 2, m))))
    mix_old_mu <- as.normmix(ps[i-1, ], mus[[i-1]], sigmas[[i-1]])
    approx_old_mu <- approx_normmix(mix_old_mu, win)
    for (j in 1:m) {
      mus[[i]][j, ] <- sample_mu(j, old_mu = mus[[i-1]][j, ],
                                 approx_old_mu[j],
                                 sigmas[[i-1]][, , j])

      sigmas[[i]][, , j] <- sample_sigma(j, mu = mus[[i]][j, ],
                                         old_sigma = sigmas[[i-1]][ , , j])
      #sample indicators zij
#       sig1 <- sigmas[[i]][, , j]
#       invsig1 <- solve(sig1)
#       if (truncate == TRUE) {
#         ###add code
#       }
#       consts[i, j] <- (1/approx[j])*det(2*pi*sig1)^(-0.5)
    }

    # sample ps
    ds <- gam + summary(as.factor(zmultinom))
    ps[i, ] <- rdirichlet(1, ds)

    mix <- as.normmix(ps[i, ], mus[[i]], sigmas[[i]])
    den <- matrix(NA_real_, n, mix$m)
    for (k in 1:mix$m) {
      den[, k] <- mvtnorm::dmvnorm(pp, mix$mus[[k]], mix$sigmas[[k]])
      den[, k] <- den[, k] * mix$ps[k]
    }

    qij <- t(apply(den, 1, function(x) x / sum(x)))
    propz <- apply(qij, 1, sample, x = 1:m, size = 1, replace = T)

    ratio <- ifelse(any(summary(as.factor(propz)) < 2), 0, 1)

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
