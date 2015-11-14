## ---- echo=FALSE---------------------------------------------------------
library(knitr)
opts_chunk$set(fig.width = 5, fig.height = 4)
library(rgl)
knit_hooks$set(webgl = hook_webgl)


## ------------------------------------------------------------------------
library(sppmix)
data(redwood)
Window(redwood)

## ---- results='hide'-----------------------------------------------------
fit <- est_mix_damcmc(redwood, 4, L = 10000) 

## ------------------------------------------------------------------------
post_mix <- get_post(fit, burnin = 2000)
post_mix

## ------------------------------------------------------------------------
plot_contour(mix = post_mix$post_normmix, lambda = post_mix$mean_lambda,  
             pattern = redwood, win = Window(redwood))

## ------------------------------------------------------------------------
summary(fit)

## ------------------------------------------------------------------------
plot_ind(fit)

## ---- webgl=TRUE---------------------------------------------------------
plot(post_mix$post_normmix, post_mix$mean_lambda, Window(fit$data))

