#' @export
FixLS_da<- function(gens, win, truemix=NULL, maxz=1, m=2, burnin=1000,
                 xlab1="x",ylab1="y", plot_result = FALSE)
{
  xlims1 <- c(win$xrange)
  ylims1 <- c(win$yrange)
  L <- dim(gens$genps)[1]
   permgens <- PostGenGetBestPerm_sppmix(gens$allgens_List)
   if (plot_result == TRUE) {
#    post_ps <- colMeans(permgens$permuted_ps[-(1:burnin), ])
#    mus <- apply(permgens$permuted_mus[, , -(1:burnin)], 1:2, mean)
#
#    mean_mat <- function(mats) Reduce("+", mats) / length(mats)
#    sigmas <- apply(permgens$permuted_sigmas[-(1:burnin), ], 2, mean_mat)
#
#    mean_lambda <- mean(gens$genlamdas[burnin:L])
#
#    post_mus <- post_sigmas <- vector("list", m)
#    for (i in 1:m) {
#      post_mus[[i]] <- mus[i, ]
#      post_sigmas[[i]] <- matrix(sigmas[, i], 2, 2)
#    }
#    post_normix = normmix(post_ps, post_mus, post_sigmas)

#    n = gens$data$n
#     par(mfrow=c(1,1))
#
#    titleLines <- list(
#      bquote(paste(lambda,"=",.(mean_lambda),", n=",.(n),", m=",.(m)," components")),
#         title1 = "Posterior mean (permutated labels)"
#    )
#    # Now output each line The text in the list is converted
#    # to expressions do.call
#
#    plot.default(gens$data,pch=20,
#         xlab=xlab1,
#         ylab=ylab1,
#         xlim=xlims1,
#         ylim=ylims1,main="")
#    mtext(do.call(expression, titleLines),side=3,line=0:1)
#
#    for(i in 1:m)
#    {
#      center = c(post_normix$mus[[i]])
#      points(center[1], center[2] ,pch=20,col="red")
#    }
#
#    plot.normmix(post_normix, mean_lambda, win = win,
#                 title1 = "Posterior mean intensity surface (permutated labels)")

#   if(!is.null(truemix))
#     ShowStats(permgens$permuted_ps,permgens$permuted_mus,truemix)
#   ShowChains(permgens$permuted_ps,permgens$permuted_mus)
#   allpermgens=list(genps = permgens$permuted_ps,
#                    genlamdas=allgens$genlamdas ,
#                    allgens = permgens$permuted_gens)
#   if(Get_User_Input_sppmix(
#     "Show average of intensity surfaces \n(slow operation, permuted realizations)?"))
#     plot_avgsurf(gens, win, burnin = burnin)
   }
   perum_fit <- list(allgens_list = permgens$permuted_gens,
                     genps = permgens$permuted_ps,
                     genmus = permgens$permuted_mus,
                     gensigmas = permgens$permuted_sigmas,
                     genzs = permgens$best_perm,
                     genlamdas = gens$genlamdas,
                     data = gens$data
                      )
   class(perum_fit) <- "damcmc_res"
   return(invisible(perum_fit))
}
