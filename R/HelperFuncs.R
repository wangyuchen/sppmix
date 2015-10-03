meansmix<- function(truemix,n=50,m=1,xlims=c(0,10),ylims=c(0,10),L=10000,burnin=1000,LL=31,trunc=FALSE)
{
#  L=5000
#  m=1
#  burnin=1000
#  n=50
#  xlims<-c(-10,10)
#  ylims<-c(-10,10)
#  mu=c(0,0)
 # sig=rbind(c(5,0),c(0,2))
  #rnorm_cpp(n,mu,sig)
 # m <- diag(2);
  #matrix(rnorm(4), 2, 2)
  #getEigenValues(m)
  #mvrnormArma(10,c(0,0),diag(2))
  # microbenchmark(DAMCMC2d(n,mu,sig))
#build a list of lists
#  x=rMultinomial_sppmix(1,rDirichlet_sppmix(rep(1,m)))

  #  truemix=GenNormalMixture(n,m,xlims,ylims,r=2,trunc)

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
  cat("\n posterior mean p returned\n",gens$meanps[i])
  cat("\n trueps\n",truemix[[i]]$p)
#  cat("\n posterior meanmus[",i,"]\n",meanmus[i,1],meanmus[i,2])
  cat("\n posterior meanmus returned[ ",i,"]\n",gens$meanmus[i,])
  cat("\n truemean[",i,"]\n",truemix[[i]]$mu)
#  cat("\n posterior mean sigma[",i,"]\n",meansigmas[,,i])
  cat("\n posterior mean sigma returned[",i,"]\n",gens$meansigmas[,,i])
  cat("\n true sigma[",i,"]\n",truemix[[i]]$sigma,"\n")
}
#  print(GetStats_sppmix(gens$genps[,1],alpha=0.05))
  
#cat("passed")
graphics.off()
while (rgl::rgl.cur()>0)
{
  rgl::rgl.close()
}

mix=vector("list", m)
for(i in 1:m)
{
  mix[[i]]=list(p = gens$meanps[i], 
             mu = gens$meanmus[i,],
             sigma = gens$meansigmas[,,i])
}
z<-PlotNormalMixture(mix1=mix,data=data,m=m,xlims=xlims,ylims=ylims,L=100)

#library(graphics)
#windows()
#persp(x = gens$ticsx, y = gens$ticsy, 
#        z = gens$AvgofPostIntensity)
#title2=paste("Mixture with",m,"components",
#    "\n",n,"points")
LL=length(gens$ticsx);
Plot3d_sppmix(x = as.vector(gens$ticsx),
             y = as.vector(gens$ticsy), 
             z = as.matrix(gens$AvgofPostIntensity))

z1=ApproxAvgPostIntensity(gens$allgens,
                          LL,burnin,
                         as.vector(gens$ticsx),
                         as.vector(gens$ticsy))

Plot3d_sppmix(x = as.vector(gens$ticsx),
              y = as.vector(gens$ticsy), 
              z = as.matrix(z1))

#cat("passedPlotNormalMixture")
windows()
par(mfrow=c(1,2))
plot(gens$genmus[1,2,],type="l")
plot(gens$genmus[1,1,],type="l")

}

GenNormalMixture<- function(n=50,m=1,  
                            xlims=c(-3,3),ylims=c(-3,3),r=20,trunc=FALSE)
{
  #  truemix=GenNormalMixture(n,m,xlims=c(0,10),ylims=c(0,10),r=5,false)
  ps<<-rDirichlet_sppmix(rep(1,m))
  mix=vector("list", m)
  for(i in 1:m)
  {
    mui=rbind(runif(1,xlims[1],xlims[2]),
              runif(1,ylims[1],ylims[2]))
    sigi=rWishart_sppmix(r,diag(2))
    compi=list(p = ps[i], mu = mui, sigma = sigi)
    mix[[i]]=compi
  }
  #  data
  #  data=rnorm2_sppmix(n,mu,sig)
  mix1<<-mix
  genmix<<-genNormMix_sppmix(n,mix1)
  data<<-genmix[[1]];
  comps<<-genmix[[2]];
  print(table(comps))
#  opar <- par()      # make a copy of current settings
  par(mfrow=c(1,1))
  #  par(pch="20") 
  #  windows()
  
  plot(data,pch=20,xlab="longitude",
       ylab="latitude",xlim=xlims, ylim=ylims,
       main=paste("Mixture with",m,"components",
                  "\n",n,"points"))
  #       ,sub=paste(n,"points"))
  for(i in 1:m)
  {
    points(mix1[[i]]$mu[1],mix1[[i]]$mu[2],pch=20,col="red")
  }
#  par(opar)          
  return(mix1)
}

PlotNormalMixture<- function(mix1,data,
                             m=1,  
                             xlims=c(-3,3),
                             ylims=c(-3,3),L=100)
{
  #  PlotNormalMixture(truemix,data,m,xlims=c(0,10),ylims=c(0,10),L=100)
  n=nrow(data) 
#  opar <- par()      # make a copy of current settings
  par(mfrow=c(1,1))
  #  par(pch="20") 
  #  windows()
  
  plot(data,pch=20,xlab="longitude",
       ylab="latitude",xlim=xlims, ylim=ylims,
       main=paste("Mixture with",m,"components",
                  "\n",n,"points"));
  #       ,sub=paste(n,"points"))
  for(i in 1:m)
  {
    points(mix1[[i]]$mu[1],mix1[[i]]$mu[2],pch=20,col="red");
  }
  
  xcoord <- seq(xlims[1], xlims[2], length.out = L);
  ycoord <- seq(ylims[1], ylims[2], length.out = L);
  
  z <- dNormMix_sppmix(mix1, xcoord, ycoord);
  
#  cat("passed")
#  cat(mix1[[1]]$p)
  # assign colors to heights for each point
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                   "#7FFF7F", "yellow", "#FF7F00", "red",
                                   "#7F0000"))
  col <- jet.colors(100)[findInterval(z, seq(min(z), max(z), length = 100))];
  
  rgl::open3d();
  rgl::persp3d(x = xcoord, y = ycoord, z = z,
               color = col, alpha = .8);
#  shapelist3d(cube3d(), x,y,z, col=grey.colors(5),add=T)
  
#  library(fields)
  ## add color bar
#  image.plot(legend.only=T, zlim=range(z), col=col)
  #  par(opar)  
  
  return(z)
}

Plot3d_sppmix<- function(xcoord,ycoord,z,title1="Poisson Process Intensity")
{
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                   "#7FFF7F", "yellow", "#FF7F00", "red",
                                   "#7F0000"))
  col <- jet.colors(100)[findInterval(z, seq(min(z), max(z), length = 100))]
  
#  title2=as.character(title1)
#  library(rgl)
#  rgl::rgl.bringtotop()
  rgl::open3d()
  rgl::persp3d(x = xcoord, y = ycoord, z = z,
               color = col, alpha = .8, xlab="x",
               ticktype="detailed",ylab="y",
               zlab="Intensity")
  rgl::title3d(main=title1)
#  rgl::axes3d(edges=c('x','y+','z'),pos=c(0,0,0),box=FALSE)
#  rgl::axis3d(edge='x',pos=c(0,0,0))
#  rgl::axis3d(edge='y',pos=c(0,0,0))
#  rgl::axis3d(edge='z',pos=c(0,0,0))
  #  rgl::decorate3d(main=title1)
}