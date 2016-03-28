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
} else {
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
  mix_demo <- rnormmix(m, sig0, sigdf, xlim = xlims, ylim = ylims)
  demo_surf <- to_int_surf(mix_demo, lambda = lamda, win = win)
  truemix <- rsppmix(demo_surf, truncate = truncated)
  plot(truemix, mus = mix_demo$mus)

if(Get_User_Input_sppmix("Simulate from the posterior (DAMCMC)?"))
  gens <- est_mix_damcmc(truemix, m = m, L = L, truncate = truncated)
  post_demo <- get_post(gens, burnin = burnin)
#  cat("passed")
if(Get_User_Input_sppmix("Show basic 2d and 3d plots?"))
{
  zmax_truemix <- max(dnormmix(demo_surf, xlim = xlims, ylim = ylims))
  zmax_genmeanmix <- max(dnormmix(post_demo, xlim = xlims, ylim = ylims))
  zmax <- max(c(zmax_truemix, zmax_genmeanmix))

  print(plotmix_2d(post_demo, pattern = truemix))
  plot(post_demo,
       zlims = c(0, 1.1*zmax), truncate = truncated,
       main = paste("Intensity surface of posterior means,",m,
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
  test_labswitch(gens)

if(Get_User_Input_sppmix("Show average of intensity surfaces (slow operation)?"))
  plot_avgsurf(gens, win = win, burnin = burnin, LL = LL, zlims=c(0,1.1*zmax))

if(Get_User_Input_sppmix("Apply relabeling algorithm?"))
{
  perm_post <- FixLS_da(gens, plot_result = TRUE)
}
#  sppmix::BDMCMC2d_sppmix(20,gendata,xlims,ylims,L,31,FALSE,1,20,c(15,.01,3,2,1,1))

if(Get_User_Input_sppmix("Run the Birth-Death MCMC fit?"))
{
  gensBD <- est_mix_bdmcmc(pp = truemix, m = maxnumcomp, truncate = truncated,
                           lambda1 = 1, lambda2 = 10,
                           hyper = c(5,.01,3,2,1,1), L = L)
if(Get_User_Input_sppmix("Show Birth-Death MCMC plots?")) {
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
    plot(postBD, zlims = c(0, 1.1*zmax), truncate = truncated,
         main = paste("Intensity surface of posterior means, MAP m = ",
                        which.max(tab),
                        "components,",truemix$n,"points"))
}
  if(Get_User_Input_sppmix("Show average surfaces of Birth-Death MCMC?")) {
    plot(gensBD)
  }
}

