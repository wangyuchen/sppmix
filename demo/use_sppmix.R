library(MASS)
# mu=c(0,0)
# sig=rbind(c(1,0),c(0,1))
#rnorm_cpp(n,mu,sig)
#  m <- diag(2)
#matrix(rnorm(4), 2, 2)
#getEigenValues(m)
#mvrnormArma(10,c(0,0),diag(2))
# microbenchmark(DAMCMC2d(n,mu,sig))
#  data
if(Get_User_Input_sppmix("Use default parameter values?"))
{
  xlims<<-c(0,10)
  ylims<<-c(0,10)
  L<<-5000
  m<<-3
  burnin<<-1000
  lamda<<-100
  r<<-30
  LL<<-50
  maxnumcomp<<-10
}else
{
  m<<-Get_User_Input_sppmix("the true number of components",modeYN=0)
  xlims<<-c(
    Get_User_Input_sppmix("minx",modeYN=0),
    Get_User_Input_sppmix("maxx",modeYN=0))
  ylims<<-c(
    Get_User_Input_sppmix("miny",modeYN=0),
    Get_User_Input_sppmix("maxy",modeYN=0))
  L<<-Get_User_Input_sppmix("number of iterations",modeYN=0)
  burnin<<-Get_User_Input_sppmix("burnin",modeYN=0)
  lamda<<-Get_User_Input_sppmix("lambda (affects number of events)",modeYN=0)
  r<<-Get_User_Input_sppmix("parameter r (affects the cov matrices of components)",modeYN=0)
  #    if(Get_User_Input_sppmix("Apply truncation to the window specified (results in slow operations)"))
  #      truncated<<-TRUE
  #    else
  #      truncated<<-FALSE
  LL<<-Get_User_Input_sppmix("grid side length (the larger the slower the operations but much prettier...)",modeYN=0)
  maxnumcomp<<-Get_User_Input_sppmix("maximum number of components for BDMCMC",modeYN=0)
}
# data=rnorm2_sppmix(n,mu,sig)
#  Plots_off()

if(Get_User_Input_sppmix("Remove all plots?"))
  Plots_off()

if(Get_User_Input_sppmix("Apply truncation?"))
  truncated<<-TRUE else  truncated<<-FALSE

if(Get_User_Input_sppmix("Generate the true mixture?"))
  truemix<<-GenNormalMixture(lamda,m,xlims,ylims,r,truncated)

if(Get_User_Input_sppmix("Simulate from the posterior (DAMCMC)?"))
  gens<<-DAMCMC2d_sppmix(gendata,xlims,ylims,m,L,LL,trunc=truncated)

#  cat("passed")
if(Get_User_Input_sppmix("Show basic 2d and 3d plots?"))
{
  zmax_truemix=lamda*GetMixtureMaxz_sppmix(truemix,
                                           100,xlims,ylims);
  mix_of_postmeans<<-#MakeMixtureList_sppmix(
    MakeMixtureList(gens$allgens_List,burnin)
  mean_lambda<<-mean(gens$genlamdas[burnin:L]);
  zmax_genmeanmix=mean_lambda *
    GetMixtureMaxz_sppmix(mix_of_postmeans,
                          100,xlims,ylims);
  #find the highest z
  maxz_height<<-max(c(zmax_truemix,zmax_genmeanmix))
  #do all the plotting with a common maximum z value
  PlotNormalMixture(mix1=truemix,data1=gendata,
                    m1=m,lamda1=lamda,xlims1=xlims,
                    ylims1=ylims,L1=100,
                    title1="True mixture",zlims1=c(0,1.1*maxz_height),
                    title3d=paste("True mixture intensity surface,",m,"components,",n,"points"))

  PlotNormalMixture(mix1=mix_of_postmeans,
                    data1=gendata,
                    m1=m,lamda1=mean_lambda,xlims1=xlims,
                    ylims1=ylims,
                    L1=100,title1="Posterior means",zlims1=c(0,1.1*maxz_height),
                    title3d=paste("Intensity surface of posterior means,",m,"components,",n,"points"))
}
if(Get_User_Input_sppmix("Show Chains and Stats?"))
{
  #  plot_ind(gens)
  ShowChains(gens$genps,gens$genmus)
  ShowStats(gens$genps,gens$genmus,truemix)
}

if(Get_User_Input_sppmix("Check for label switching?"))
  CheckLabels(gens$genmus[,,burnin:L])

if(Get_User_Input_sppmix("Show average of intensity surfaces (slow operation)?"))
  Show3dAvgofsurfaces(gens,LL,burnin,xlims,ylims,zlims=c(0,1.1*maxz_height))

if(Get_User_Input_sppmix("Apply relabeling algorithm?"))
{
  FixLabels(gens,data1=gendata,truemix,1.1*maxz_height,
            m1=m,xlims1=xlims,ylims1=ylims)
}
#  sppmix::BDMCMC2d_sppmix(20,gendata,xlims,ylims,L,31,FALSE,1,20,c(15,.01,3,2,1,1))

if(Get_User_Input_sppmix("Run the Birth-Death MCMC fit?"))
{
  gensBD<<-BDMCMC2d_sppmix(maxnumcomp,gendata,xlims,ylims,L,LL,FALSE,1,10,c(5,.01,3,2,1,1))
  if(Get_User_Input_sppmix("Show Birth-Death MCMC plots?"))
    PostGenBDMCMC_sppmix(gensBD,maxz_height)
}

