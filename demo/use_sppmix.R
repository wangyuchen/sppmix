library(MASS)
if(Get_User_Input_sppmix("Use default parameter values?"))
{
  xlims <- c(0,10)
  ylims <- c(0,10)
  win <- square(10)
  L <- 5000
  m <- 3
  burnin <- 1000
  lamda <- 100
  sig0 <- 0.1
  sigdf <- 5
  LL <- 50
  maxnumcomp <- 10
}else
{
  m <- Get_User_Input_sppmix("the true number of components",modeYN=0)
  xlims <- c(
    Get_User_Input_sppmix("minx",modeYN=0),
    Get_User_Input_sppmix("maxx",modeYN=0))
  ylims <- c(
    Get_User_Input_sppmix("miny",modeYN=0),
    Get_User_Input_sppmix("maxy",modeYN=0))
  win <- owin(xlims,ylims)
  L <- Get_User_Input_sppmix("number of iterations",modeYN=0)
  burnin <- Get_User_Input_sppmix("burnin",modeYN=0)
  lamda <- Get_User_Input_sppmix("lambda (affects number of events)",modeYN=0)
  sig0 <- Get_User_Input_sppmix(
    "Tunning parameter in generating random matrix from Wishart distribution)",
    modeYN=0)
  sigdf <- Get_User_Input_sppmix(
    "Degree of freedom in generating random matrix from Wishart distribution.",
    modeYN = 0)
  #    if(Get_User_Input_sppmix("Apply truncation to the window specified (results in slow operations)"))
  #      truncated <- TRUE
  #    else
  #      truncated <- FALSE
  LL <- Get_User_Input_sppmix("grid side length (the larger the slower the operations but much prettier...)",modeYN=0)
  maxnumcomp <- Get_User_Input_sppmix("maximum number of components for BDMCMC",modeYN=0)
}
# data=rnorm2_sppmix(n,mu,sig)
#  Plots_off()

if(Get_User_Input_sppmix("Remove all plots?"))
  Plots_off()

if(Get_User_Input_sppmix("Apply truncation?"))
  truncated <- TRUE else  truncated <- FALSE

if(Get_User_Input_sppmix("Generate the true mixture?"))
  #truemix <- GenNormalMixture(lamda,m,xlims,ylims,r,truncated)
  mix_demo <- rnormmix(m, sig0, sigdf, win)
  truemix <- rsppmix(lamda, mix = mix_demo, truncate = truncated, win = win)
  true_mus <- data.frame(matrix(unlist(mix_demo$mus), m, 2, byrow = TRUE))
  names(true_mus) <- c("x","y")
  plot(truemix)+ggplot2::geom_point(data = true_mus, size = 3, colour ="red")

if(Get_User_Input_sppmix("Simulate from the posterior (DAMCMC)?"))
  gens <- est_mix_damcmc(truemix, m = m, L = L, truncate = truncated)
  post_demo <- get_post(gens, burnin = burnin)
#  cat("passed")
if(Get_User_Input_sppmix("Show basic 2d and 3d plots?"))
{
  zmax_truemix <- max(lamda*dnormmix(mix_demo,win));
  zmax_genmeanmix <- max(post_demo$mean_lambda*dnormmix(post_demo$post_normmix,
                                                        win))
  zmax <- max(c(zmax_truemix, zmax_genmeanmix))

  plot(mix = mix_demo, lambda = lamda, win = win,
       zlims = c(0, 1.1*zmax), truncate = truncated,
       title1 = paste("Intensity surface of the true model with,",m,
                      "components,",truemix$n,"points"))
  plot(post_demo$post_normmix, post_demo$mean_lambda, win = win,
       zlims = c(0, 1.1*zmax), truncate = truncated,
       title1 = paste("Intensity surface of posterior means,",m,
                      "components,",truemix$n,"points"))
#   PlotNormalMixture(mix1=mix_of_postmeans,
#                     data1=gendata,
#                     m1=m,lamda1=mean_lambda,xlims1=xlims,
#                     ylims1=ylims,
#                     L1=100,title1="Posterior means",zlims1=c(0,1.1*maxz_height),
#                     title3d=paste("Intensity surface of posterior means,",m,
#                                   "components,", truemix$n,"points"))
}
if(Get_User_Input_sppmix("Show Chains and Stats?"))
{
  #  plot_ind(gens)
  plot_chains(gens, burnin = burnin)
  summary(gens)
}

if(Get_User_Input_sppmix("Check for label switching?"))
  test_labswitch(gens$genmus[,,burnin:L])

if(Get_User_Input_sppmix("Show average of intensity surfaces (slow operation)?"))
  plot_avgsurf(gens, win = win, burnin = burnin, LL = 60, zlims=c(0,1.1*zmax))

if(Get_User_Input_sppmix("Apply relabeling algorithm?"))
{
  perm_post <- FixLS_da(gens, plot_result = TRUE)
}
#  sppmix::BDMCMC2d_sppmix(20,gendata,xlims,ylims,L,31,FALSE,1,20,c(15,.01,3,2,1,1))

if(Get_User_Input_sppmix("Run the Birth-Death MCMC fit?"))
{
  gensBD <- est_mix_bdmcmc(pp = truemix, m = m, truncate = truncated,
                           lambda = 1, lambdab = 10,
                           hyper = c(5,.01,3,2,1,1), L = L)
  if(Get_User_Input_sppmix("Show Birth-Death MCMC plots?"))
    cat("Frequency table for the number of components")
  print(table(gensBD$numcomp))
  tab=tabulate(gensBD$numcomp,nbins=gensBD$maxnumcomp)
    mp <- barplot(tab,names.arg=1:gensBD$maxnumcomp,
                 xlab="Number of Components",ylab="Iterations",
                 main="Distribution of the number of components"
                 ,ylim=c(0,1.2*max(tab)))
    plot(gensBD$numcomp, xlab="Iteration", ylab="Number of components",
         type="l", main="Generated chain for the number of components")
    postBD <- get_post(gensBD, comp =which.max(tab), burnin = burnin)
    plot(postBD$post_normmix, postBD$mean_lambda, win = win,
         zlims = c(0, 1.1*zmax), truncate = truncated,
         title1 = paste("Intensity surface of posterior means, MAP m=,",m,
                        "components,",truemix$n,"points"))
}

