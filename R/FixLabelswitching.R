#' Fix label switching
#'
#' Fix label switching by apply the best permutation to the posterior
#' realization (MCMC result)
#' @inheritParams plot_avgsurf
#' @param xlab1 The label for x-axis
#' @param ylab1 The label for y-axis
#' @param plot_result Logical flag indicating whether plot the point pattern
#' and intensity surface after permutation. The default is FALSE.
#' @author Athanasios Christou Micheas, Jiaxun Chen, Yuchen Wang
#' @export
#' @examples
#' # generate data
#' mix2 <- normmix(ps=c(.4, .6), mus=list(c(0.3, 0.3), c(0.7, 0.7)),
#' sigmas = list(.02*diag(2), .01*diag(2)))
#' pp2 <- rsppmix(100,mix2,square(1))
#' # Run Data augmentation MCMC and get posterior realizations
#' post = est_mix_damcmc(pp2,L = 5000,2,truncate = F)
#' # get posterior mean for each parameter
#' post_mean = get_post(post)
#' # plot the estimated intensity surface
#' plot(post_mean$post_normmix, post_mean$mean_lambda, square(1))
#' # Fix label switching
#' post_fixed = FixLS_da(post, plot_result = TRUE)
FixLS_da<- function(fit, burnin = length(fit$allgens_List) / 10,
                 xlab1 = "x",ylab1 = "y", plot_result = FALSE)
{
  win <- domain(fit$data)
  m <- dim(fit$genmus)[1]
  xlims1 <- c(win$xrange)
  ylims1 <- c(win$yrange)
  L <- dim(fit$genps)[1]
   permgens <- PostGenGetBestPerm_sppmix(fit$allgens_List)
   if (plot_result == TRUE) {
   post_ps <- colMeans(permgens$permuted_ps[-(1:burnin), ])
   mus <- apply(permgens$permuted_mus[, , -(1:burnin)], 1:2, mean)

   mean_mat <- function(mats) Reduce("+", mats) / length(mats)
   sigmas <- apply(permgens$permuted_sigmas[-(1:burnin), ], 2, mean_mat)

   mean_lambda <- mean(gens$genlamdas[burnin:L])

   post_mus <- post_sigmas <- vector("list", m)
   for (i in 1:m) {
     post_mus[[i]] <- mus[i, ]
     post_sigmas[[i]] <- matrix(sigmas[, i], 2, 2)
   }
   post_normix = normmix(post_ps, post_mus, post_sigmas)

   n = fit$data$n
    par(mfrow=c(1,1))

   titleLines <- list(
     bquote(paste(lambda,"=",.(mean_lambda),", n=",.(n),", m=",.(m)," components")),
        title1 = "Posterior mean (permutated labels)"
   )
   # Now output each line The text in the list is converted
   # to expressions do.call

   plot.default(fit$data,pch=20,
        xlab=xlab1,
        ylab=ylab1,
        xlim=xlims1,
        ylim=ylims1,main="")
   mtext(do.call(expression, titleLines),side=3,line=0:1)

   for(i in 1:m)
   {
     center = c(post_normix$mus[[i]])
     points(center[1], center[2] ,pch=20,col="red")
   }

   plot.normmix(post_normix, mean_lambda, win = win,
                title1 = "Posterior mean intensity surface (permutated labels)")

#   if(!is.null(truemix))
#     ShowStats(permgens$permuted_ps,permgens$permuted_mus,truemix)
#   ShowChains(permgens$permuted_ps,permgens$permuted_mus)
#   allpermgens=list(genps = permgens$permuted_ps,
#                    genlamdas=allgens$genlamdas ,
#                    allgens = permgens$permuted_gens)
#   if(Get_User_Input_sppmix(
#     "Show average of intensity surfaces \n(slow operation, permuted realizations)?"))
#     plot_avgsurf(fit, win, burnin = burnin)
   }
   perum_fit <- list(allgens_List = permgens$permuted_gens,
                     genps = permgens$permuted_ps,
                     genmus = permgens$permuted_mus,
                     gensigmas = permgens$permuted_sigmas,
                     genlamdas = fit$genlamdas,
                     data = fit$data
                      )
   class(perum_fit) <- "damcmc_res"
   return(invisible(perum_fit))
}

#' Test if posterior realizations of mu have label switching
#'
#' Test if there is a label switching of the posterior realizations of mu by
#' testing if the mean of each component changed dramatically during MCMC.
#'
#' @param genmus posterior realizations from DAMCMC algorithm.
#' @author Athanasios (Sakis) Micheas
#' @export
test_labswitch<- function(genmus) {
  m <- dim(genmus)[1]
  cat("\nChecking for label switching...\n")
  for (i in 1:m) {
    if(Check4LabelSwitching_sppmix(genmus[i, 1, ])) {
      cat("Label switching present. \nPermute the labels to get a better fit,
          \nor obtain the average of the surfaces\n ")
      return(TRUE)
    }
    if(Check4LabelSwitching_sppmix(genmus[i, 2, ])) {
      cat("Label switching present. \nPermute the labels to get a better fit,
          \nor obtain the average of the surfaces\n ")
      return(TRUE)
    }
  }
  cat("No Label switching detected")
  return(invisible(FALSE))
}

