# rbivnormal <- function(n, mu = 0, Sigma = diag(2)) {
#   z <- rnorm(2 * n)
#   decomp <- chol(Sigma)
#   x <- mu + t(decomp) %*% matrix(z, 2, n)
#   return(t(x))
# }
#
# dbivnormal <- function(x, mu = 0, Sigma = diag(2)) {
#   exp(-.5 * t(x - mu) %*% solve(Sigma) %*% (x - mu)) / sqrt(2 * pi * det(Sigma))
# }
#
# rwishart <- function(n, Sigma) {
#   x <- rbinorm(n, 0, Sigma)
#   xxp <- list()
#   for (i in 1:n) {
#     xxp[[i]] <- x[i, ] %*% t(x[i, ])
#   }
#   return(Reduce("+", xxp))
# }

rdirichlet <- function (n, alpha) {
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  return(x/as.vector(sm))
}
