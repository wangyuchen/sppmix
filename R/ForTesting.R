#package creation summary
#Do Clean+Rebuild, test, then check, then B+R
#' @export
meansmix<- function(data,truemix,lamda=50,m=1,xlims=c(0,10),ylims=c(0,10),L=10000,burnin=1000,LL=31,trunc=FALSE)
{
#  L=5000
#  m=1
#  burnin=1000
#  lamda=50
#  xlims<-c(-10,10)
#  ylims<-c(-10,10)
#  mu=c(0,0)
 # sig=rbind(c(5,0),c(0,2))
  #rnorm_cpp(n,mu,sig)
 # m <- diag(2);
  #matrix(rnorm(4), 2, 2)
  #getEigenValues(m)
  #mvrnormArma(10,c(0,0),diag(2))
  #library(microbenchmark)
# microbenchmark(sppmix::rUniform_sppmix(5000),sppmix::rUniform_sppmix1(5000),runif(5000),times=1000)

#build a list of lists
#  x=rMultinomial_sppmix(1,rDirichlet_sppmix(rep(1,m)))

  #truemix=GenNormalMixture(n,m,xlims,ylims,r,trunc)

  gens<<-DAMCMC2d_sppmix(data,xlims,ylims,m,L,burnin,LL,trunc)
 #gens$allgens[[10]][[2]]#10nth realization, 2nd mix comp
if(0)#no need, the returned means work
{
meanmus=matrix(0,m,2)
meansigmas=array(0,c(2,2,m))
# meaninvsigmas=zeros(m,2,2);
meanp=matrix(0,m,1)
for (i in 1:m)
{
  meanp[i]=mean(gens$genps[(burnin+1):L,i]);
  meanmus[i,1]=mean(gens$genmus[i,1,(burnin+1):L]);
  meanmus[i,2]=mean(gens$genmus[i,2,(burnin+1):L]);
  sigs=matrix(unlist(gens$gensigmas[(burnin+1):L,i]),L-burnin,4,byrow=TRUE)
  meansigmas[1,1,i]=mean(sigs[,1]);
  meansigmas[1,2,i]=mean(sigs[,2]);
  meansigmas[2,1,i]=mean(sigs[,3]);
  meansigmas[2,2,i]=mean(sigs[,4]);
}
}
# x[1,1,1]
#  mix1 <- normmix(ps=c(.3, .7), mus=list(c(0, 0), c(1, 1)),sigmas=list(.01*diag(2), .01*diag(2)))
#cat("\n data avg \n",mean(data[,1]),mean(data[,2]))
#windows()
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
  poststats=GetStats_sppmix(gens$genps[,i],alpha=0.05)
  cat("\n----------------Component ",i,"------------\n")
  cat(paste("Probability: true =",truemix[[i]]$p))
  cat(paste("\nProbability: posterior mean =",poststats$Mean,"\n"))
  cat(paste(poststats$CredibleSetConfidence, "% Credible Set:\n[",poststats$CredibleSet[1],
            ",",poststats$CredibleSet[2],"]" ))
  poststats=GetStats_sppmix(gens$genmus[i,1,],alpha=0.05)
  cat(paste("\nMean vector, x-coord: true =",truemix[[i]]$mu[1]))
  cat(paste("\nMean vector, x-coord: post mean =",poststats$Mean,"\n"))
  cat(paste(poststats$CredibleSetConfidence, "% Credible Set:\n[",poststats$CredibleSet[1],
            ",",poststats$CredibleSet[2],"]" ))
  poststats=GetStats_sppmix(gens$genmus[i,2,],alpha=0.05)
  cat(paste("\nMean vector, x-coord: true =",truemix[[i]]$mu[2]))
  cat(paste("\nMean vector, x-coord: post mean =",poststats$Mean,"\n"))
  cat(paste(poststats$CredibleSetConfidence, "% Credible Set:\n[",poststats$CredibleSet[1],
            ",",poststats$CredibleSet[2],"]" ))

}
cat("\n----------------Component stats done------------\n")
cat("\nNOTE: if you have the truth, then what is called
component 1 might be called component 2 in the fit.
This is not label switching.\n")

Plots_off()

#cat("passed")

PlotNormalMixture(mix1=truemix,data=data,
                     m=m,lamda=lamda,xlims=xlims,
                     ylims=ylims,L=100,
                     title1=paste(
            "True mixture intensity,",
                m,"components,",n,"points"))

mix=MakeMixtureList(gens);

z<-PlotNormalMixture(mix1=mix,data=data,
   m=m,lamda=gens$meanlamda,xlims=xlims,ylims=ylims,
   L=100,title1=paste("Mixture intensity of posterior means,",m,"components,",n,"points"))

#library(graphics)
#windows()
#persp(x = gens$ticsx, y = gens$ticsy,
#        z = gens$AvgofPostIntensity)
#title2=paste("Mixture with",m,"components",
#    "\n",n,"points")
LL=length(gens$ticsx);

ApproxAvgPostIntensityz<<-ApproxAvgPostIntensity(
  gens$allgens,gens$genlamdas,LL,burnin,
  as.vector(gens$ticsx),as.vector(gens$ticsy))

Plot3d_sppmix(x = as.vector(gens$ticsx),
              y = as.vector(gens$ticsy),
              z = ApproxAvgPostIntensityz,
              title1=paste("Average of",L-burnin,"posterior realizations of the intensity function")
              ,xlims,ylims,zlims)

#Plot3d_sppmix(x = as.vector(gens$ticsx),
#              y = as.vector(gens$ticsy),
#              z = z1,
#              title1="Posterior fit using ApproxAvgPostIntensity")

#cat("passedPlotNormalMixture")
if(1)
{
  windows()
  par(mfrow=c(m,1))
  plot(gens$genps[,1],xlab="Iteration",ylab="p",
    type="l",main=
    "Generated mixture Probabilities\nComponent 1")
  for (i in 2:m)
  {
    plot(gens$genps[,i],xlab="Iteration",ylab="p",
         type="l",main=
           paste("Generated mixture Probabilities\nComponent",i))
  }

  windows()
  par(mfrow=c(m,2))
  plot(gens$genmus[1,1,],xlab="Iteration",ylab=
         bquote(mu),
       type="l",main=
         "Generated mixture means\nComponent 1, x-coord")
  plot(gens$genmus[1,2,],xlab="Iteration",ylab=
         bquote(mu),type="l",main=
         "Generated mixture means\nComponent 1, y-coord")
  for (i in 2:m)
  {
    plot(gens$genmus[i,1,],xlab="Iteration",ylab=
           bquote(mu),type="l",main=
           paste("Generated mixture means\nComponent"
                 ,i,", x-coord"))
    plot(gens$genmus[i,2,],xlab="Iteration",ylab=
           bquote(mu),type="l",main=
           paste("Generated mixture means\nComponent"
                 ,i,", y-coord"))
  }
}


#cat(gens$allgens[[L]][[1]]$p)

}

#' @export
GenNormalMixture<- function(lamda=50,m=1,
    xlims=c(-3,3),ylims=c(-3,3),r=20,trunc=FALSE,
    xlab1="x",ylab1="y")
{
  #  truemix=GenNormalMixture(lamda,m,xlims=c(0,10),ylims=c(0,10),r=5,false)
  X11()
  ps<<-rDirichlet_sppmix(rep(1,m))
  mix=vector("list", m)
  Rangex=xlims[2]-xlims[1];
  Rangey=ylims[2]-ylims[1];
  Amat=rbind(c(1/Rangex^2,0),c(0,1/Rangey^2))
  for(i in 1:m)
  {
    mui=rbind(rUnifab_sppmix(xlims[1],xlims[2]),
              rUnifab_sppmix(ylims[1],ylims[2]))
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

  titleLines <- list(
    bquote(paste(lambda,"=",.(lamda),", n=",.(n))),
    bquote(paste("True Mixture has ",.(m)," components"))
  )
#  Now output each line The text in the list is converted to expressions do.call

  plot(gendata,pch=20,
       xlab=xlab1,#xlab="longitude",
       ylab=ylab1,#ylab="latitude",
       xlim=xlims,
       ylim=ylims,main="")
  mtext(do.call(expression, titleLines),side=3,line=0:1)
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
PlotNormalMixture<- function(mix1,data,
                             m=1,lamda=100,
                             xlims=c(-3,3),
                             ylims=c(-3,3),L=100,
                             title1="Mixture intensity",
                             xlab1="x",ylab1="y",
                             zlims=c(0,1))
{
  #  z=sppmix::PlotNormalMixture(truemix,gendata,m,lamda,xlims=c(0,10),ylims=c(0,10),L=100,title1="Mixture",zlims=c(0,1))
  windows()
  n=nrow(data)
#  opar <- par()      # make a copy of current settings
  par(mfrow=c(1,1))
  #  par(pch="20")
  #  windows()

  titleLines <- list(
    bquote(paste(lambda,"=",.(lamda),", n=",.(n),", m=",.(m)," components")),
#    bquote(paste("Mixture Intensity with ",.(m)," components"))
    title1
  )
  #  Now output each line The text in the list is converted to expressions do.call

  plot(data,pch=20,
       xlab=xlab1,#xlab="longitude",
       ylab=ylab1,#ylab="latitude",
       xlim=xlims,
       ylim=ylims,main="")
  mtext(do.call(expression, titleLines),side=3,line=0:1)
  #  bquote("Mixture Intensity, m =" ~.(m)~","
  #        ~lambda ==.(lamda)~", n " == .(n)))

  for(i in 1:m)
  {
    points(mix1[[i]]$mu[1],mix1[[i]]$mu[2],pch=20,col="red")
  }

#  plot(data,pch=20,xlab="longitude",
#       ylab="latitude",xlim=xlims, ylim=ylims,
#       main=paste("Mixture with",m,"components",
 #                 "\n",n,"points"));

  xcoord <- seq(xlims[1], xlims[2], length.out = L);
  ycoord <- seq(ylims[1], ylims[2], length.out = L);

  zcoord <- lamda*dNormMix_sppmix(mix1, xcoord, ycoord);

  Plot3d_sppmix(xcoord,ycoord,zcoord,
                title1=paste(title1,"intensity"),
                xlims,ylims,zlims)
  #,pos1=c(xlims[1],ylims[2],0))
  #              paste("Mixture with",
   #            m,"components,",n,"points"))
  return(zcoord)
}

#' @export
Plot3d_sppmix<- function(xcoord,ycoord,z,
      title1="Poisson Process Intensity",
      xlims=c(0,1),ylims=c(0,1),zlims=c(0,1))
{
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                   "#7FFF7F", "yellow", "#FF7F00", "red",
                                   "#7F0000"))
  col <- jet.colors(100)[findInterval(z, seq(min(z), max(z), length = 100))]

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
#              ,userMatrix =
#              rgl::rotationMatrix(pi/4,1, 0, 0))
  U=rgl::par3d("userMatrix")
  rgl::par3d(userMatrix=
            rgl::rotate3d(U,pi/4,0,0,1))
  zmax=max(z)
  Rangez=zmax-min(z);
  rgl::persp3d(x = xcoord, y = ycoord, z = z,
     color = col, xlab="x",ylab="y",zlab="",
     zlim=c(zlims[1]-0.01*Rangez,zlims[2]+0.01*Rangez),
     box = FALSE, axes = FALSE)
  rgl::axis3d('x')
  rgl::axis3d('y')
  rgl::axis3d('z-+'#-=low x-coord, +=high y-coord
#              ,zlim=c(zlims[1]-0.01*Rangez,zlims[2]+0.01*Rangez)
              ,pos = c(xlims[1], ylims[2], 0))
#  Rangex=max(xcoord)-min(xcoord);
#  Rangey=max(ycoord)-min(ycoord);
#  Rangez=max(z)-min(z);
#  zmax=max(z)
#  rgl::observer3d(-xmax, -ymax, zmax,auto = TRUE)
  rgl::title3d(main=NULL)
  rgl::text3d(xlims[2],ylims[2],
              zmax+0.2*Rangez,texts=title1)
#  rgl::text3d(xmax+pos[1],ymax+pos[2],
#              zmax+pos[3],texts=title1)
#  rgl::axes3d(edges=c('x','y+','z'),pos=c(0,0,0),box=FALSE)
#  rgl::axis3d(edge='x',pos=c(0,0,0))
#  rgl::axis3d(edge='y',pos=c(0,0,0))
#  rgl::axis3d(edge='z',pos=c(0,0,0))
  #  rgl::decorate3d(main=title1)
}

#' @export
Plots_off<- function()
{
graphics.off()
while (rgl::rgl.cur()>0)
{
  rgl::rgl.close()
}
}

#' @export
Go<- function()
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
xlims<<-c(0,10)
ylims<<-c(0,10)
L<<-5000
m<<-3
burnin<<-1000
lamda<<-100
r<<-30
truncated<<-FALSE
LL<<-50
# data=rnorm2_sppmix(n,mu,sig)
#  Plots_off()
 truemix<<-GenNormalMixture(lamda,m,xlims,ylims,r,truncated)
#  PlotNormalMixture(truemix,gendata,m,lamda,
#                    xlims,ylims,L=100)
 gens<<-DAMCMC2d_sppmix(gendata,xlims,ylims,m,L,burnin,LL,trunc=FALSE)

  Plots_off()

#  if(0){

#  cat("passed")

  zmax_truemix=lamda*GetMixtureMaxz_sppmix(truemix);
  mix=MakeMixtureList(gens)
  zmax_genmeanmix=gens$meanlamda*GetMixtureMaxz_sppmix(mix);

  maxz_height<<-max(c(zmax_truemix,zmax_genmeanmix))
  PlotNormalMixture(mix1=truemix,data=gendata,
      m=m,lamda=lamda,xlims=xlims,ylims=ylims,L=100,
      title1="True mixture",zlims=c(0,maxz_height))
                      #paste(
                      #"True Mixture Intensity,",
                      #m,"components,",n,"points"))

  PlotNormalMixture(mix1=mix,data=gendata,
     m=m,lamda=gens$meanlamda,xlims=xlims,ylims=ylims,
     L=100,title1="Posterior means",
     zlims=c(0,maxz_height))
                         #paste("Mixture Intensity of posterior means,",m,"components,",n,"points"))

#  Show3dPlots(gens)
if(0){
  plot_ind(gens)
  ShowChains(gens)
  ShowStats(gens,truemix)
  Show3dPlots(gens)
  FixLabels(gens,truemix)
}
#  sppmix::BDMCMC2d_sppmix(20,gendata,xlims,ylims,L,burnin,31,FALSE,1,20,c(15,.01,3,2,1,1))
if(1)
{
  maxnumcomp<<-10
  gensBD<<-BDMCMC2d_sppmix(maxnumcomp,gendata,xlims,ylims,L,burnin,LL,FALSE,1,10,c(5,.01,3,2,1,1))
  PostGenBDMCMC_sppmix(gensBD)
}
  #  setwd("D:/sppmix/images")

}


#' @export
PostGenBDMCMC_sppmix<- function(gensBD)
{#processes output from a Birth-Death fit
  cat("Frequency table for the number of components")
  tab=table(gensBD$numcomp)
  print(tab)

#  probk=tab
#  meank=sum(probk.*[1:1:maxnumcomp]);
#  [val1,mapk]=max(probk);
  windows()
  hist(gensBD$numcomp,breaks=0:gensBD$maxnumcomp,xlab="Component",labels=TRUE,main="Distribution of the number of components",xlim=c(0.5,gensBD$maxnumcomp),ylim=c(0,1.2*max(table(gensBD$numcomp))))
  windows()
  plot(gensBD$numcomp,xlab="Iteration",ylab="Number of components",type="l",main="Generated chain for the number of components")
  gridvals=GetGrid_sppmix(LL,c(xlims[1],ylims[1]),c(xlims[2],ylims[2]));
  ticsx<<-as.vector(gridvals[[1]]);
  ticsy<<-as.vector(gridvals[[2]]);
  distr_numcomp<<-GetCompDistr_sppmix(gensBD$numcomp[burnin:L],maxnumcomp)
  BayesianModelAvgIntensity<<-ApproxBayesianModelAvgIntensity_sppmix(
    gensBD$allgens_List[burnin:L],
    gensBD$genlamdas[burnin:L],
    gensBD$numcomp[burnin:L],
    distr_numcomp,
    0,maxnumcomp,
#    1,5,
    ticsx,ticsy)
  Plot3d_sppmix(xcoord = ticsx,ycoord = ticsy,
                z = BayesianModelAvgIntensity,
                title1=paste("Bayesian model average of",L-burnin,"posterior realizations")
                ,xlims,ylims,zlims=c(0,maxz_height))

}

#' @export
Show3dPlots<- function(gens)#,LL=50,xlims=c(0,10),ylims=c(0,10))
{
#  sppmix::Show3dPlots(gens,LL,xlims,ylims)
  #library(graphics)
  #windows()
  #persp(x = gens$ticsx, y = gens$ticsy,
  #        z = gens$AvgofPostIntensity)
  #title2=paste("Mixture with",m,"components",
  #    "\n",n,"points")
  LL=length(gens$ticsx);
  #Plot3d_sppmix(x = as.vector(gens$ticsx),
  #             y = as.vector(gens$ticsy),
  #             z = gens$AvgofPostIntensity,
  #             title1=paste("Average of",L-burnin,"posterior realizations of the intensity function"))

#  gridvals=GetGrid_sppmix(LL,c(xlims[1],ylims[1]),c(xlims[2],ylims[2]));
#  ticsx=gridvals[[1]];
#  ticsy=gridvals[[2]];
  ApproxAvgPostIntensityz<<-ApproxAvgPostIntensity(
    gens$allgens,gens$genlamdas,LL,burnin,
    as.vector(gens$ticsx),as.vector(gens$ticsy))
#  ApproxAvgPostIntensityz<<-ApproxAvgPostIntensity(
#    gens$allgens,gens$genlamdas,LL,burnin,
#    as.vector(ticsx),as.vector(ticsy));

  Plot3d_sppmix(xcoord = as.vector(gens$ticsx),
                ycoord = as.vector(gens$ticsy),
                z = ApproxAvgPostIntensityz,
                title1=paste("Average of",L-burnin,"posterior realizations of the intensity function")
                ,xlims,ylims,zlims)

 # ShowChains(gens)

}

#' @export
ShowChains<- function(gens)
{
    windows()
    par(mfrow=c(m,1))
    plot(gens$genps[,1],xlab="Iteration",ylab="p",
         type="l",main=
           "Generated mixture probabilities\nComponent 1")
    for (i in 2:m)
    {
      plot(gens$genps[,i],xlab="Iteration",ylab="p",
           type="l",main=
             paste("Generated mixture probabilities\nComponent",i))
    }

    windows()
    par(mfrow=c(m,2))
    plot(gens$genmus[1,1,],xlab="Iteration",ylab=
           bquote(mu),
         type="l",main=
           "Generated mixture means\nComponent 1, x-coord")
    plot(gens$genmus[1,2,],xlab="Iteration",ylab=
           bquote(mu),type="l",main=
           "Generated mixture means\nComponent 1, y-coord")
    for (i in 2:m)
    {
      plot(gens$genmus[i,1,],xlab="Iteration",ylab=
             bquote(mu),type="l",main=
             paste("Generated mixture means\nComponent"
                   ,i,", x-coord"))
      plot(gens$genmus[i,2,],xlab="Iteration",ylab=
             bquote(mu),type="l",main=
             paste("Generated mixture means\nComponent"
                   ,i,", y-coord"))
    }

}

#' @export
FixLabels<- function(allgens,truemix,burnin=1000)
{
#  z=sppmix::FixLabels(gens,truemix)
#  ShowChains(allgens)
  permgens<<-PostGenGetBestPerm_sppmix(allgens$allgens);
  allpermgens=allgens;
  allpermgens[[1]]=permgens$permuted_gens;
  allpermgens[[2]]=permgens$permuted_ps;
  allpermgens[[3]]=permgens$permuted_mus;
  allpermgens[[4]]=permgens$permuted_sigmas;
  permuted_means=GetAllMeans_sppmix(permgens$permuted_gens,burnin)
  allpermgens[[6]]=permuted_means$meanps;
  allpermgens[[7]]=permuted_means$meanmus;
  allpermgens[[8]]=permuted_means$meansigmas;
  mix=MakeMixtureList(allpermgens);
  z<-PlotNormalMixture(mix1=mix,data=gendata,
                       m=m,lamda=gens$meanlamda,xlims=xlims,ylims=ylims,
                       L=100,title1="Posterior mean (permutated labels)")

  ShowStats(allpermgens,truemix)
  ShowChains(allpermgens)
  Show3dPlots(allpermgens)
}

#' @export
MakeMixtureList<- function(gens)
{
  #takes the DAMCMC output and makes a mixture list
  #based on its means, sppmix::MakeMixtureList(gens)
  m=length(gens$meanps);
  mix=vector("list", m);
  for(i in 1:m)
  {
    mix[[i]]=list(p = gens$meanps[i],
                mu = gens$meanmus[i,],
                sigma = gens$meansigmas[,,i]);
  }
  return (mix)
}

#' @export
ShowStats<- function(gens,truemix)
{
  #sppmix::ShowStats(gens,truemix)
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
  poststats=GetStats_sppmix(gens$genps[,i],alpha=0.05)
  cat("\n----------------Component ",i,"------------\n")
  cat(paste("Probability: true =",truemix[[i]]$p))
  cat(paste("\nProbability: posterior mean =",poststats$Mean,"\n"))
  cat(paste(poststats$CredibleSetConfidence, "% Credible Set:\n[",poststats$CredibleSet[1],
            ",",poststats$CredibleSet[2],"]" ))
  poststas=GetStats_sppmix(gens$genmus[i,1,],alpha=0.05)
  cat(paste("\nMean vector, x-coord: true =",truemix[[i]]$mu[1]))
  cat(paste("\nMean vector, x-coord: post mean =",poststats$Mean,"\n"))
  cat(paste(poststats$CredibleSetConfidence, "% Credible Set:\n[",poststats$CredibleSet[1],
            ",",poststats$CredibleSet[2],"]" ))
  poststats=GetStats_sppmix(gens$genmus[i,2,],alpha=0.05)
  cat(paste("\nMean vector, x-coord: true =",truemix[[i]]$mu[2]))
  cat(paste("\nMean vector, x-coord: post mean =",poststats$Mean,"\n"))
  cat(paste(poststats$CredibleSetConfidence, "% Credible Set:\n[",poststats$CredibleSet[1],
            ",",poststats$CredibleSet[2],"]" ))

}
cat("\n----------------Component stats done------------\n")
#cat("\nNOTE: if you have the truth, then what is called component 1 might be called component 2 in the fit. This is not label switching. Check the plots of the generated chains for erratic behavior (sudden jumps) and if present, run the FixLabels routine to find the best permutation.\n")
}
