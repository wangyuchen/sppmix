## ---- message=FALSE, warning=FALSE---------------------------------------
library(sppmix)
int_surf <- normmix(ps = c(.3, .2, .5),
                    mus = list(c(0.2, 0.2), c(0.3, 0.5), c(.8, .8)),
                    sigmas = list(.01*diag(2), .01*diag(2), .01*diag(2)),
                    lambda = 100,
                    win = square(1))

## ----webgl=TRUE----------------------------------------------------------
plot(int_surf)

## ------------------------------------------------------------------------
pp <- rsppmix(int_surf)
plot(pp, int_surf$mus)
fit <- est_mix_damcmc(pp, m=3, truncate = TRUE)

## ------------------------------------------------------------------------
postmean <- get_post(fit)
plot(postmean)

