#package creation summary
#Do Clean+Rebuild, test, then check, then B+R
#Written by Sakis Micheas, 2015

#' @export
GenNormalMixture<- function(lamda=50,m=1,
    xlims=c(-3,3),ylims=c(-3,3),r=20,trunc=FALSE,
    xlab1="x",ylab1="y")
{
  #  truemix=GenNormalMixture(lamda,m,xlims=c(0,10),ylims=c(0,10),r=5,false)
#  windows()
  ps<<-rDirichlet_sppmix(rep(1,m))
  mix=vector("list", m)
  Rangex=xlims[2]-xlims[1];
  Rangey=ylims[2]-ylims[1];
  Amat=rbind(c(1/Rangex^2,0),c(0,1/Rangey^2))
  for(i in 1:m)
  {
    mui=rbind(runif(1,xlims[1],xlims[2]),
              runif(1,ylims[1],ylims[2]))
    sigi=rWishart_sppmix(r,Amat)
    compi=list(p = ps[i], mu = mui, sigma = sigi)
    mix[[i]]=compi
  }
  #  data
  #  data=rnorm2_sppmix(n,mu,sig)
  mix1<-mix
  genmix<<-rNormMix_sppmix(lamda,mix1)
  gendata<<-genmix[[1]];
  n<<-nrow(gendata);
  truecomps<<-genmix[[2]];
  print(table(truecomps))
#  opar <- par()      # make a copy of current settings
  par(mfrow=c(1,1))
  #  par(pch="20")
  #  windows()

#  library(grDevices)
#  require(graphics)
if(0)
{
  count=n;
  if(trunc)
  {
    count=0;
    indices=rep(0,n);
    for(i in 1:n)
    {
      if(gendata[i,1]>=xlims[1]
         && gendata[i,1]<=xlims[2]
         && gendata[i,2]>=ylims[1]
         && gendata[i,2]<=ylims[2])
      count=count+1;
      indices[i]=1;
    }
    cat(paste(n-count,'points were outside W=[',
      xlims[1],',',xlims[2],']x[',ylims[1],',',
      ylims[2],']'))
#    data=data[(indices==1),];
  }
}
  checkin=CheckInWindow_sppmix(gendata,xlims,ylims,trunc);
  count=checkin$count_inW;

  titleLines <- list(
    bquote(paste("Displaying ",.(count)," points, ",.(n-count)," points outside window")),
    bquote(paste(lambda,"=",.(lamda),", n=",.(n))),
    bquote(paste("True Mixture has ",.(m)," components"))
  )
#  Now output each line The text in the list is converted to expressions do.call

  plot(gendata,pch=20,
       xlab=xlab1,#xlab="longitude",
       ylab=ylab1,#ylab="latitude",
       xlim=xlims,
       ylim=ylims,main="")
  mtext(do.call(expression, titleLines),side=3,line=0:2)
#  bquote("Mixture Intensity, m =" ~.(m)~","
#        ~lambda ==.(lamda)~", n " == .(n)))

  for(i in 1:m)
  {
    points(mix1[[i]]$mu[1],mix1[[i]]$mu[2],pch=20,col="red")
  }
#  par(opar)
  return(mix1)
}

#' @export
PlotNormalMixture<- function(mix1,data1,
                             m1=1,lamda1=100,
                             xlims1=c(-3,3),
                             ylims1=c(-3,3),L1=100,
                             title1="Mixture intensity",
                             xlab1="x",ylab1="y",
                             zlims1=c(0,1),
                             title3d="")
{
 # windows()
  n=nrow(data1)
#  opar <- par()      # make a copy of current settings
  par(mfrow=c(1,1))
  #  par(pch="20")
  #  windows()

  titleLines <- list(
    bquote(paste(lambda,"=",.(lamda1),", n=",.(n),", m=",.(m1)," components")),
#    bquote(paste("Mixture Intensity with ",.(m)," components"))
    title1
  )
  #  Now output each line The text in the list is converted to expressions do.call

  plot(data1,pch=20,
       xlab=xlab1,#xlab="longitude",
       ylab=ylab1,#ylab="latitude",
       xlim=xlims1,
       ylim=ylims1,main="")
  mtext(do.call(expression, titleLines),side=3,line=0:1)
  #  bquote("Mixture Intensity, m =" ~.(m)~","
  #        ~lambda ==.(lamda)~", n " == .(n)))

  for(i in 1:m1)
  {
    points(mix1[[i]]$mu[1],mix1[[i]]$mu[2],pch=20,col="red")
  }

#  plot(data,pch=20,xlab="longitude",
#       ylab="latitude",xlim=xlims, ylim=ylims,
#       main=paste("Mixture with",m,"components",
 #                 "\n",n,"points"));

  xcoord <- seq(xlims1[1], xlims1[2], length.out = L1);
  ycoord <- seq(ylims1[1], ylims1[2], length.out = L1);

  if(0==1)
  {#use the package routine
    normmix1=MakeNormMixFromMixtureList(mix1);
    plot(normmix1, lamda,
        win=spatstat::owin(xlims1, ylims1),
        L = 100,
      title1="Poisson Process Intensity", truncate = TRUE)
  }
  else
  {#use the Plot3d_sppmix routine
    zcoord <- lamda*dNormMix_sppmix(mix1, xcoord, ycoord);
    Plot3d_sppmix(xcoord,ycoord,zcoord,
                  title1=title3d,
                  xlims1,ylims1,zlims1)
  }
  #,pos1=c(xlims[1],ylims[2],0))
  #              paste("Mixture with",
   #            m,"components,",n,"points"))
#  return(zcoord)
}


#' @export
Plot2d_sppmix<- function(data1,lamda1=100,
                 xlims1=c(-3,3),
                 ylims1=c(-3,3),
                 title1="Mixture intensity",
                 xlab1="x",ylab1="y")
{
  # windows()
  n=nrow(data1)
  #  opar <- par()      # make a copy of current settings
  par(mfrow=c(1,1))
  #  par(pch="20")
  #  windows()

  titleLines <- list(
    bquote(paste(lambda,"=",.(lamda1),", n=",.(n))),
    #    bquote(paste("Mixture Intensity with ",.(m)," components"))
    title1
  )
  #  Now output each line The text in the list is converted to expressions do.call

  plot(data1,pch=20,
       xlab=xlab1,#xlab="longitude",
       ylab=ylab1,#ylab="latitude",
       xlim=xlims1,
       ylim=ylims1,main="")
  mtext(do.call(expression, titleLines),side=3,line=0:1)
}

#' @export
Plot3d_sppmix<- function(xcoord,ycoord,zcoord,
      title1="Poisson Process Intensity",
      xlims=c(0,1),ylims=c(0,1),zlims=NULL)
{
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                   "#7FFF7F", "yellow", "#FF7F00", "red",
                                   "#7F0000"))
  cols <- jet.colors(100)[findInterval(zcoord, seq(min(zcoord), max(zcoord), length = 100))]

  if(is.null(zlims))
    zlims=c(0,max(zcoord))
#  title2=as.character(title1)
#  library(rgl)
#  rgl::rgl.bringtotop()
#  rgl::par3d(params=rgl::r3dDefaults)
  #.Platform
  height=300;
  width=500;
  if(.Platform$OS.type=='windows')
  {
    scr_width <- system("wmic desktopmonitor get screenwidth", intern=TRUE)
    scr_height <- system("wmic desktopmonitor get screenheight", intern=TRUE)
    height=as.numeric(scr_height[length(scr_height)-1])
    width=as.numeric(scr_width[length(scr_width)-1])
  }

  #  height=as.integer(system("wmic desktopmonitor get screenheight"));
#  width=as.integer(system("wmic desktopmonitor get screenwidth"));
#  cat(width/2.0)
  rgl::open3d(windowRect=c(width/5,
                           height/7,
                           4*width/5,
                           6*height/7),
              zoom=1.2)
#rotation about the x-axis, 45 degrees
  U=rgl::par3d("userMatrix")
  rgl::par3d(userMatrix=
            rgl::rotate3d(U,pi/4,0,0,1))
  zmax=max(zcoord)
  Rangez=zmax-min(zcoord);
  rgl::persp3d(x = xcoord, y = ycoord, z = zcoord,
        color = cols, xlab="x",ylab="y",zlab="",
        zlim=c(zlims[1]-0.01,zlims[2]),
        box = FALSE, axes = FALSE)
  rgl::axis3d('x')
  rgl::axis3d('y')
  rgl::axis3d('z-+',pos = c(xlims[1], ylims[2], 0))
  rgl::rgl.lines(c(xlims[1], xlims[1]),
                 c(ylims[2], ylims[2]), zlims, color = 'black')
  rgl::title3d(main=NULL)
  rgl::text3d(xlims[2],ylims[2],
              zlims[2]#+0.2*Rangez
              ,texts=title1)
  rgl::bgplot3d(suppressWarnings(
    fields::image.plot(legend.only = TRUE,
                       smallplot= c(.8,.82,0.05,.7),
                       zlim = zlims,
                       col = jet.colors(100))))
}

#' @export
Plot3dGrayScale_sppmix<- function(xcoord,ycoord,zcoord,
                                  title1="Poisson Process Intensity",
                                  xlims=c(0,1),ylims=c(0,1),zlims=NULL)
{
  cols <- gray.colors(100)[findInterval(zcoord, seq(min(zcoord), max(zcoord), length = 100))]
  height=300;
  width=500;
  if(.Platform$OS.type=='windows')
  {
    scr_width <- system("wmic desktopmonitor get screenwidth", intern=TRUE)
    scr_height <- system("wmic desktopmonitor get screenheight", intern=TRUE)
    height=as.numeric(scr_height[length(scr_height)-1])
    width=as.numeric(scr_width[length(scr_width)-1])
  }
  rgl::open3d(windowRect=c(width/5,
                           height/7,
                           4*width/5,
                           6*height/7),
              zoom=1.2)
  U=rgl::par3d("userMatrix")
  rgl::par3d(userMatrix=
               rgl::rotate3d(U,pi/4,0,0,1))
  zmax=max(zcoord)
  Rangez=zmax-min(zcoord);
  rgl::persp3d(x = xcoord, y = ycoord, z = zcoord,
               color = cols, xlab="x",ylab="y",zlab="",
               zlim=c(zlims[1]-0.01,zlims[2]),
               box = FALSE, axes = FALSE)
  rgl::axis3d('x')
  rgl::axis3d('y')
  rgl::axis3d('z-+',pos = c(xlims[1], ylims[2], 0))
  rgl::rgl.lines(c(xlims[1], xlims[1]),
                 c(ylims[2], ylims[2]), zlims, color = 'black')
  rgl::title3d(main=NULL)
  rgl::text3d(xlims[2],ylims[2],
              zlims[2]#+0.2*Rangez
              ,texts=title1)
  rgl::bgplot3d(suppressWarnings(
    fields::image.plot(legend.only = TRUE,
      smallplot= c(.8,.82,0.05,.7),
      zlim = zlims,col = gray.colors(100))))
}


#' @export
Save_AllOpenRglGraphs<- function(
  dir1="D:/sppmix/images",filename1="RglGraph")
{
  setwd(dir1)
#  graphics.off()
  count=1;
  while (rgl::rgl.cur()>0)
  {
    rgl::rgl.snapshot(paste(filename1,count,".png",sep = ""))
    rgl::rgl.close()
    count=count+1;
  }
  if (count==2)
    cat(paste("Saved 1 plot in",dir1))
  else
    cat(paste("Saved",count-1,"plots in",dir1))
}

#' @export
Go<- function()
{
#just for testing, run Demo_sppmix instead
  #do all the plotting with a common maximum z value
  xlims<<-c(0,10)
  ylims<<-c(0,10)
  L<<-5000
  m<<-5
  burnin<<-1000
  lamda<<-100
  r<<-30
  truncated<<-FALSE
  LL<<-50
  maxnumcomp<<-10

  truemix<<-GenNormalMixture(lamda,m,xlims,ylims,r,truncated)

  gens<<-DAMCMC2d_sppmix(gendata,xlims,ylims,m,L,trunc=FALSE)

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

  PlotNormalMixture(mix1=truemix,data1=gendata,
                    m1=m,lamda1=lamda,
                    xlims1=xlims,ylims1=ylims,L1=100,
                    title1="True mixture",
                    zlims1=c(0,1.1*maxz_height),
                    title3d=paste("True mixture intensity surface,",m,"components,",n,"points"))

  PlotNormalMixture(mix1=mix_of_postmeans,
                    data1=gendata,
                    m1=m,lamda1=mean_lambda,
                    xlims1=xlims,ylims1=ylims,
                    L1=100,
                    title1="Posterior means",
                    zlims1=c(0,1.1*maxz_height),
                    title3d=paste("Intensity surface of posterior means,",m,"components,",n,"points"))

  gensBD<<-BDMCMC2d_sppmix(maxnumcomp,gendata,xlims,ylims,L,LL,FALSE,1,10,c(5,.01,3,2,1,1))
  PostGenBDMCMC_sppmix(gensBD,maxz_height)

}

#' @export
Demo_sppmix<- function()
{
  library(microbenchmark)
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
    m<<-1
    burnin<<-1000
    lamda<<-100
    r<<-30
    LL<<-50
    maxnumcomp<<-10
  }
  else
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
    truncated<<-TRUE
  else
    truncated<<-FALSE

  if(Get_User_Input_sppmix("Generate the true mixture?"))
    truemix<<-GenNormalMixture(lamda,m,xlims,ylims,r,truncated)

  if(Get_User_Input_sppmix("Simulate from the posterior (DAMCMC)?"))
    gens<<-DAMCMC2d_sppmix(gendata,xlims,ylims,m,L,trunc=truncated)

#  class(gens) <<- "damcmc_res"
 #   cat("passed")
  if(Get_User_Input_sppmix("Show basic 2d and 3d plots?"))
  {
    if (m>1)
      plot_ind(gens)
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
    ShowChains(gens$genps,gens$genmus,m=m)
#    cat("passed")
    ShowStats(gens$genps,gens$genmus,truemix)
  }

  if(m>1 && Get_User_Input_sppmix("Check for label switching?"))
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

}


#' @export
PostGenBDMCMC_sppmix<- function(gensBD,maxz_height)
{#processes output from a Birth-Death fit
  cat("Frequency table for the number of components")
  print(table(gensBD$numcomp))
  tab=tabulate(gensBD$numcomp,nbins=gensBD$maxnumcomp)

#  probk=tab
#  meank=sum(probk.*[1:1:maxnumcomp]);
#  [val1,mapk]=max(probk);
#  windows()
#  windows()
  mp<<-barplot(tab,names.arg=1:gensBD$maxnumcomp,
          xlab="Number of Components",ylab="Iterations",
          main="Distribution of the number of components"
          ,ylim=c(0,1.2*max(tab)))
  for(i in 1:gensBD$maxnumcomp)
    text(x=mp,y=tab+0.1*max(tab),labels=tab)
  #mtext(side=1,at=c(mp[i,1],-1),text=paste(tab[i]))
  #  mtext(side = 1, at = colMeans(mp), line = 1,
 #       text = paste("Mean"), col = "red")

  #  windows()
  plot(gensBD$numcomp,xlab="Iteration",ylab="Number of components",type="l",main="Generated chain for the number of components")
  distr_numcomp<<-GetCompDistr_sppmix(gensBD$numcomp[burnin:L],maxnumcomp)

  MAPcompList=GetMax_sppmix(tab)
  MAPcomp=MAPcompList$pos+1;
  mapgens<<-GetBDCompRealiz_sppmix(
    gensBD$allgens_List[burnin:L],
    gensBD$genlamdas[burnin:L],
    gensBD$numcomp[burnin:L],MAPcomp)

  mix_of_postmeans<<-MakeMixtureList_sppmix(
    mapgens$newgenBD,0)
  mean_lambda=mean(mapgens$newlamdas);
#  maxz_height=mean_lambda*
#    GetMixtureMaxz_sppmix(mix_of_postmeans,
#                          100,xlims,ylims);
#  cat(maxz_height)
#  cat("passed")
  PlotNormalMixture(mix1=mix_of_postmeans,
    data1=gendata,m1=MAPcomp,lamda1=mean_lambda,
    xlims1=xlims,ylims1=ylims,L1=100,
    title1="Posterior means",
    zlims1=c(0,1.1*maxz_height),
    title3d=paste("Intensity surface of posterior means, MAP m=",m,",",n,"points"))

  if(0)
  {
    permgens<<-PostGenGetBestPerm_sppmix(mapgens$newgenBD);
    #  permuted_means=GetAllMeans_sppmix(permgens$permuted_gens,burnin)
    mix_of_permuted_means<<-MakeMixtureList_sppmix(permgens$permuted_gens,0);
    PlotNormalMixture(mix1=mix_of_permuted_means,
      data1=gendata,m1=MAPcomp,lamda1=mean_lambda,
      xlims1=xlims,ylims1=ylims,L1=100,
      title1="Posterior mean (permutated labels)",
      zlims1=c(0,1.1*maxz_height),
      title3d = "Posterior mean intensity surface (permutated labels)")
  }
  if(1)
  {

  if(Get_User_Input_sppmix("Compute Bayesian model average (slow operation)?"))
  {
    BayesianModelAvgIntensity<<-
      ApproxBayesianModelAvgIntensity_sppmix(
      gensBD$allgens_List[burnin:L],
      gensBD$genlamdas[burnin:L],
      gensBD$numcomp[burnin:L],
      distr_numcomp,1,maxnumcomp,LL,xlims,ylims)
    gridvals=GetGrid_sppmix(LL,xlims,ylims);
    ticsx<<-as.vector(gridvals[[1]]);
    ticsy<<-as.vector(gridvals[[2]]);
    Plot3d_sppmix(xcoord = ticsx,ycoord = ticsy,
                zcoord = BayesianModelAvgIntensity,
                title1=paste("Bayesian model average of",L-burnin,"posterior realizations")
                ,xlims,ylims,zlims=c(0,1.1*maxz_height))
  }
}
}

#' @export
Show3dAvgofsurfaces<- function(gens,LL=30,burnin=1000,xlims=c(0,10),ylims=c(0,10),zlims=c(0,1))
{
  L=dim(gens$genps)[1]
  #library(graphics)
  #windows()
  #persp(x = gens$ticsx, y = gens$ticsy,
  #        z = gens$AvgofPostIntensity)
  #title2=paste("Mixture with",m,"components",
  #    "\n",n,"points")
#  LL=length(gens$ticsx);
  #Plot3d_sppmix(x = as.vector(gens$ticsx),
  #             y = as.vector(gens$ticsy),
  #             z = gens$AvgofPostIntensity,
  #             title1=paste("Average of",L-burnin,"posterior realizations of the intensity function"))

  gridvals=GetGrid_sppmix(LL,xlims,ylims);
  ticsx=gridvals[[1]];
  ticsy=gridvals[[2]];
#  ApproxAvgPostIntensityz<<-ApproxAvgPostIntensity(
#    gens$allgens,gens$genlamdas,LL,burnin,
#    as.vector(gens$ticsx),as.vector(gens$ticsy))
  ApproxAvgPostIntensityz<<-ApproxAvgPostIntensity(
    gens$allgens,gens$genlamdas,LL,burnin,
    xlims,ylims);

  Plot3d_sppmix(xcoord = as.vector(ticsx),
                ycoord = as.vector(ticsy),
                zcoord = ApproxAvgPostIntensityz,
                title1=paste("Average of",L-burnin,"posterior realizations of the intensity function")
                ,xlims,ylims,zlims)


}

#' @export
ShowChains<- function(genps,genmus,m=5)
{
 #   windows()
#    par(mfrow=c(min(3,m),1))
    plot(genps[,1],xlab="Iteration",ylab="p",
         type="l",main=
           "Generated mixture probabilities\nComponent 1")
    if(m>1)
    for (i in 2:min(3,m))
    {
      plot(genps[,i],xlab="Iteration",ylab="p",
           type="l",main=
             paste("Generated mixture probabilities\nComponent",i))
    }

#    windows()
#    par(mfrow=c(min(3,m),2))
    plot(genmus[1,1,],xlab="Iteration",ylab=
           bquote(mu),
         type="l",main=
           "Generated mixture means\nComponent 1, x-coord")
    plot(genmus[1,2,],xlab="Iteration",ylab=
           bquote(mu),type="l",main=
           "Generated mixture means\nComponent 1, y-coord")

    if(m>1)
    for (i in 2:min(3,m))
    {
      plot(genmus[i,1,],xlab="Iteration",ylab=
             bquote(mu),type="l",main=
             paste("Generated mixture means\nComponent"
                   ,i,", x-coord"))
      plot(genmus[i,2,],xlab="Iteration",ylab=
             bquote(mu),type="l",main=
             paste("Generated mixture means\nComponent"
                   ,i,", y-coord"))
    }

}

#' @export
FixLabels<- function(allgens,data1,truemix=NULL,maxz=1,m1=5,xlims1=c(0,10),ylims1=c(0,10),burnin=1000)
{
#  z=sppmix::FixLabels(gens,truemix)
  permgens<<-PostGenGetBestPerm_sppmix(allgens$allgens_List);
#  permuted_means=GetAllMeans_sppmix(permgens$permuted_gens,burnin)
  mix_of_permuted_means=MakeMixtureList(permgens$permuted_gens,burnin);
  mean_lambda<<-mean(allgens$genlamdas[burnin:L]);
  PlotNormalMixture(mix1=mix_of_permuted_means,
    data1=data1,m1=m1,lamda1=mean_lambda,
    xlims1=xlims1,ylims1=ylims1,
    L1=100,title1="Posterior mean (permutated labels)",
    zlims1=c(0,maxz),title3d = "Posterior mean intensity surface (permutated labels)")

  if(!is.null(truemix))
   ShowStats(permgens$permuted_ps,permgens$permuted_mus,truemix)
  ShowChains(permgens$permuted_ps,permgens$permuted_mus,m=m)
  allpermgens=list(genps = permgens$permuted_ps,
       genlamdas=allgens$genlamdas ,
       allgens = permgens$permuted_gens)
  if(Get_User_Input_sppmix("Show average of intensity surfaces \n(slow operation, permuted realizations)?"))
    Show3dAvgofsurfaces(allpermgens,LL,burnin,xlims1,ylims1,zlims=c(0,maxz))

  return(permgens)
}

#' @export
CheckLabels<- function(genmus)
{
  cat("\nChecking for label switching...\n")
  for (i in 1:m)
  {
    if(Check4LabelSwitching_sppmix(genmus[i,1,]))
    {
      cat("Label switching present. \nPermute the labels to get a better fit,\nor obtain the average of the surfaces")
      return(TRUE);
    }
    if(Check4LabelSwitching_sppmix(genmus[i,2,]))
    {
      cat("Label switching present. \nPermute the labels to get a better fit,\nor obtain the average of the surfaces")
      return(TRUE);
    }
  }
  cat("No Label switching detected")
  return(FALSE);

}

#' @export
ShowStats<- function(genps,genmus,truemix)
{
  #sppmix::ShowStats(genps,genmus,truemix)
# windows()
#  meansmix(gendata,truemix,n,m,xlims,ylims,L,burnin,LL=51,trunc=TRUE)
#  meansmix(gendata,truemix,n,m,xlims,ylims,L,burnin,LL,trunc)
for (i in 1:m)
{
  #  cat("\n posterior meanp\n",meanp[i])
  #  cat("\n posterior mean p returned\n",gens$meanps[i])
  #  cat("\n trueps\n",truemix[[i]]$p)
  #  cat("\n posterior meanmus[",i,"]\n",meanmus[i,1],meanmus[i,2])
  #  cat("\n posterior meanmus returned[ ",i,"]\n",gens$meanmus[i,])
  #  cat("\n truemean[",i,"]\n",truemix[[i]]$mu)
  #  cat("\n posterior mean sigma[",i,"]\n",meansigmas[,,i])
  #  cat("\n posterior mean sigma returned[",i,"]\n",gens$meansigmas[,,i])
  #  cat("\n true sigma[",i,"]\n",truemix[[i]]$sigma,"\n")
  #true value and credible sets
  poststats=GetStats_sppmix(genps[,i],alpha=0.05)
  cat("\n----------------Component ",i,"------------\n")
  cat(paste("Probability: true =",truemix[[i]]$p))
  cat(paste("\nProbability: posterior mean =",poststats$Mean,"\n"))
  cat(paste(poststats$CredibleSetConfidence, "% Credible Set:\n[",poststats$CredibleSet[1],
            ",",poststats$CredibleSet[2],"]" ))
  poststas=GetStats_sppmix(genmus[i,1,],alpha=0.05)
  cat(paste("\nMean vector, x-coord: true =",truemix[[i]]$mu[1]))
  cat(paste("\nMean vector, x-coord: post mean =",poststats$Mean,"\n"))
  cat(paste(poststats$CredibleSetConfidence, "% Credible Set:\n[",poststats$CredibleSet[1],
            ",",poststats$CredibleSet[2],"]" ))
  poststats=GetStats_sppmix(genmus[i,2,],alpha=0.05)
  cat(paste("\nMean vector, x-coord: true =",truemix[[i]]$mu[2]))
  cat(paste("\nMean vector, x-coord: post mean =",poststats$Mean,"\n"))
  cat(paste(poststats$CredibleSetConfidence, "% Credible Set:\n[",poststats$CredibleSet[1],
            ",",poststats$CredibleSet[2],"]" ))

}
cat("\n----------------Component stats done------------\n")
#cat("\nNOTE: if you have the truth, then what is called component 1 might be called component 2 in the fit. This is not label switching. Check the plots of the generated chains for erratic behavior (sudden jumps) and if present, run the FixLabels routine to find the best permutation.\n")
}

#' @export
MakeNormMixFromMixtureList<- function(mix)
{
  #takes a mixture list and returns
  #a normmix object
  m=length(mix);
  ps=vector("numeric", m);
  mus=vector("list", m);
  sigmas=vector("list", m);
  for(i in 1:m)
  {
    ps[i]=as.numeric(mix[[i]]$p)
    mus[[i]]=as.vector(mix[[i]]$mu)
    sigmas[[i]]=as.matrix(mix[[i]]$sigma)
  }
  norm_mix=normmix(ps, mus, sigmas);
  return (norm_mix)
}

#' @export
MakeMixtureList<- function(allgens_List,burnin=1000)
{
  #takes the DAMCMC output and makes a mixture list
  #based on its means, sppmix::MakeMixtureList(gens)
  post_means=GetAllMeans_sppmix(allgens_List,burnin);
  m=length(post_means$meanps);
  mix=vector("list", m);
  for(i in 1:m)
  {
    mix[[i]]=list(p = post_means$meanps[i],
                  mu = post_means$meanmus[i,],
                  sigma = post_means$meansigmas[,,i]);
  }
  return (mix)
}

#' @export
Get_User_Input_sppmix<- function(prompt_string="",modeYN=1)
{
  options(warn=-1)
  #modeYN=1, ask for yes or no
  #modeYN=0, ask for a value
  if(modeYN)
  {
    ret=0;#no, 1 is yes
    while(1)
    {
      ANSWER <- readline(paste(prompt_string,
              "(Y)es or (N)o? or (Q)uit "))
      if (substr(ANSWER, 1, 1) == "n"
          || substr(ANSWER, 1, 1) == "N")
      {
        ret=0
        break
      }
      if (substr(ANSWER, 1, 1) == "y"
          || substr(ANSWER, 1, 1) == "Y")
      {
        ret=1
        break
      }
      if (substr(ANSWER, 1, 1) == "q"
          || substr(ANSWER, 1, 1) == "Q")
      {
        stop("Execution ended by the user")
      }
    }
  }
  else
  {
    ret=0;#default return value
    while(1)
    {
      #message("Enter ",prompt_string,":")
      val <- readline(paste("Enter",prompt_string,": "))
      #scan(what=double())
#      check=is.na(as.numeric(val))
      if(is.na(as.numeric(val)))
        message("Enter a number, not letters")
      else
      {
        ret=as.numeric(val)
        break
      }
    }
  }
  options(warn=0)
  return(ret)
}
