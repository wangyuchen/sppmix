## ---- echo=FALSE---------------------------------------------------------
library(knitr)
opts_chunk$set(fig.width = 6, fig.height = 4)
library(rgl)
knit_hooks$set(webgl = hook_webgl)

## ---- message=FALSE, warning=FALSE---------------------------------------
library(sppmix)
int_surf <- normmix(ps = c(.3, .3, .4),
                    mus = list(c(0.2, 0.2), c(0.3, 0.5), c(.6, .5)),
                    sigmas = list(.01*diag(2), .01*diag(2), .01*diag(2)),
                    lambda = 100,
                    win = square(1))

## ----webgl=TRUE----------------------------------------------------------
plot(int_surf, main = "Intensity surface for the true model")

## ---- message=FALSE, warning=FALSE---------------------------------------
pp <- rsppmix(int_surf)
plot(pp, int_surf$mus)
fit <- est_mix_damcmc(pp, m=3, truncate = TRUE)

## ----webgl=TRUE----------------------------------------------------------
postmean <- get_post(fit)
plot(postmean, main = "Posterior intensity surface")
plotmix_2d(postmean, fit$data)

## ------------------------------------------------------------------------
plot_chains(fit)

## ------------------------------------------------------------------------
test_labswitch(fit$genmus)

## ----message=FALSE, warning=FALSE,webgl=TRUE-----------------------------
fit2 <- FixLS_da(fit)
postmean2 <- get_post(fit2)
plot(postmean2, main = "Posterior intensity surface after fixing label switching")
plotmix_2d(postmean2, pattern = fit$data)

## ---- message=FALSE, warning=FALSE,webgl=TRUE----------------------------
plot_avgsurf(fit)

## ----message=FALSE, warning=FALSE----------------------------------------
gofts <- mc_gof(pp, postmean2, alpha = 0.05)
gofts

## ---- message=FALSE, warning=FALSE, results='hide'-----------------------
best_fit <- selectMix(pp, 1:5)

## ------------------------------------------------------------------------
best_fit

## ---- message=FALSE, warning=FALSE---------------------------------------
fitbd <- est_mix_bdmcmc(pp, 5)
tb <- tabulate(fitbd$numcomp, fitbd$maxnumcomp)
mp <- barplot(tb,names.arg=1:fitbd$maxnumcomp,
                 xlab="Number of Components",ylab="Iterations",
                 main="Distribution of the number of components"
                 ,ylim=c(0,1.2*max(tb)))

## ---- message=FALSE, warning=FALSE,webgl=TRUE----------------------------
modelbest <- get_post(fitbd, num_comp = which.max(tb))
plot(modelbest, main = paste("Instensity surface for", which.max(tb), "component model."))
plotmix_2d(modelbest, fit$data)

## ---- message=FALSE, warning=FALSE,webgl=TRUE----------------------------
plot(fitbd)

