est_mix_intensity <- function(pattern, win, m, L = 10000, burnin = 2000,
                              truncate = TRUE) {
  sample_mu <- function(j, zmultinom, pp, old_mu, sigma,
                        kappa, ksi, truncate) {
    # sample mu from proposal for one component
    sum1 <- sum(zmultinom == j)
    if (sum1 > 0) {
      newmu <- colMeans(pp[zmultinom == j,])
    } else {
      newmu <- c(0, 0)
    }

    #sample mus
    sig1 <- sigmas[[i-1]][, , j]
    invsig1 <- solve(sig1)
    cov1 <- solve(sum1*invsig1 + kappa)
    mu1 <- cov1 %*% (sum1*invsig1 %*% newmu + kappa %*% ksi)
    propmu <- mvtnorm::rmvnorm(1, mu1, cov1)
    if (truncate == TRUE) {
      ### add code
    } else {
      ratio <- 1
    }

    if (runif(1) < ratio) new_mu <- propmu else new_mu <- old_mu
    return(new_mu)
  }

  sample_sigma <- function(j, zmultinom, pp, mu, old_sigma,
                           a, beta, truncate) {
    # sample sigma from proposal distribution for one component
    sum1 <- sum(zmultinom == j)
    sumxmu <- 0
    if (sum1 > 0) {
      for (r in 1:sum1) {
        sumxmu <- sumxmu +
          t(as.matrix(pp[zmultinom == j,][r,] - mus[[i]][j,])) %*%
          (as.matrix(pp[zmultinom == j,][r,] - mus[[i]][j,]))
      }
    }

    ps2 <- solve(2*beta[ , , 1] + sumxmu)
    invsig11 <- rWishart(1, 2*a + sum1, ps2)
    propsigma <- solve(invsig11[, , 1])

    if (truncate == TRUE) {
      ### add code
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
  hmat <- a/(100*g)*kappainv
  consts <- matrix(0, L, m)

  #starting value
  mus <- list(pp[sample(1:n, size = m, replace = T),])
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
    # calculate inverse sigma
    invsigmas <- array(dim = c(2, 2, m))
    for (j in 1:m) {
      invsigmas[, , j] <- solve(sigmas[[i-1]][, , j])
    }
    sumsig <- apply(invsigmas, 1:2, sum)
    ps1 <- solve(2*hmat + 2*sumsig)
    beta <- rWishart(1, 2*g + 2*n*a, ps1)
    ds <- rep(0,  m)
    approx <- rep(1, m)

    #sample mus and sigmas
    mus <- append(mus, list(matrix(NA, m, 2)))
    sigmas <- append(sigmas, list(array(NA, dim = c(2, 2, m))))


    for (j in 1:m) {
      mus[[i]][j, ] <- sample_mu(j, zmultinom, pp,
                                 old_mu = mus[[i-1]][j, ],
                                 sigmas[[i-1]][, , j],
                                 kappa, ksi, truncate)

      sigmas[[i]][, , j] <- sample_sigma(j, zmultinom, pp,
                                         mu = mus[[i]][j, ],
                                         old_sigma = sigmas[[i-1]][ , , j],
                                         a, beta, truncate)

      ds[j] <- gam + sum(zmultinom == j)

      #sample indicators zij
      sig1 <- sigmas[[i]][, , j]
      invsig1 <- solve(sig1)
      if (truncate == TRUE) {
        ###add code
      }
      consts[i, j] <- (1/approx[j])*det(2*pi*sig1)^(-0.5)
    }

    ps[i, ]=rdirichlet(1, ds)
    mix <- as.normmix(ps[i, ],mus[[i]],sigmas[[i]])
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
    }
  }
  close(pb)
  postmus <- Reduce("+",mus[-(1:burnin)])/(L - burnin)
  postsigmas <- Reduce("+",sigmas[-(1:burnin)])/(L - burnin)
  postps <- colMeans(ps[-(1:burnin), ])
  post.mix <- as.normmix(postps,postmus,postsigmas)
  return(post.mix)
}
