## ---- message=FALSE------------------------------------------------------
library(sppmix)
ps <- c(.3, .7)
mus <- list(c(0, 0), c(1, 1))
sigmas <- list(.01 * diag(2), .01 * diag(2))
mix1 <- normmix(ps, mus, sigmas)
mix1

## ------------------------------------------------------------------------
summary(mix1)

## ------------------------------------------------------------------------
mix2 <- rnormmix(3, sig0 = .01, sigdf = 5, win = square(5))
mix2

# generate random number of components
mix3 <- rnormmix(8, sig0 = .01, sigdf = 10, square(5), rand_m = TRUE)
mix3

## ------------------------------------------------------------------------
pp1 <- rsppmix(200, mix1, square(1))
pp1

## ---- fig.width=4, fig.height=6------------------------------------------
plot(pp1)

