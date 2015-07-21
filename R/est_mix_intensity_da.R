est_mix_intensity <- function(pattern, win, m, L, burnin, truncate) {
  pattern <- cells
  win <- square(1)
  m <- 5
  L=10000
  burnin <- 2000
  truncate <- FALSE

  pp <- as.data.frame(pattern)
  n <- npoints(pattern)
  if (truncate==TRUE) {
    pattern <- pattern[spatstat::inside.owin(pattern$x,pattern$y,win)]
  }
  if (L <= burnin) {
    stop("wrong L or burnin")
  }
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
  ## start main mcmc ##
  for (i in 2:L) {
    #sample B matrix
    # calculate inverse sigma
    invsigmas <- array(dim = c(2, 2, m))
    for (j in 1:m) {
      invsigmas[, , j] <- solve(sigmas[[i-1]][, , j])
    }
    sumsig <- apply(invsigmas, 1:2, sum)
    ps1 <- solve(2*hmat + 2*sumsig)
    beta <- rWishart(1, 2*g + 2*n*a, ps1)

    #sample mus and sigmas
    mus <- append(mus, list(matrix(NA, m, 2)))
    sigmas <- append(sigmas, list(array(NA, dim = c(2, 2, m))))

    zmultinom <- sample(1:m, size = n, replace = T)
    for (j in 1:m) {
      # calculate component centers
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
      propmu <- mvtnorm::rmvnorm(1,mu1,cov1)
      if (truncate == TRUE) {
        ### add code
      } else {
        ratio <- 1
      }
      if (runif(1) < ratio) {
        mus[[i]][j,] <- propmu
      } else {
        mus[[i]][j,] <- mus[[i-1]][j,]
      }

      #sample sigmas
      if (sum1 > 0) {
        sumxmu <- 0
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

      if (runif(1) < ratio) {
        sigmas[[i]][, , j] <- propsigma
      } else {
        sigma[[i]][, , j] <- sigma[[i-1]][, , j]
      }
    }
  }
}
