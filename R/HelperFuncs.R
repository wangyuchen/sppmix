meansmix<- function(truemix,n=50,m=1,xlims=c(-3,3),ylims=c(-3,3),L=10000,burnin=1000,trunc=FALSE)
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

  gens<<-DAMCMC2d_sppmix(data,xlims,ylims,m,L,burnin,trunc)
  
meanmus=matrix(0,m,2)
meansigmas=array(0,c(m,2,2))
# meaninvsigmas=zeros(m,2,2);
meanp=matrix(0,m,1)
for (i in 1:m)
{
  meanp[i]=mean(gens$genps[(burnin+1):L,i]);
  meanmus[i,1]=mean(gens$genmus[i,1,(burnin+1):L]);
  meanmus[i,2]=mean(gens$genmus[i,2,(burnin+1):L]);
  meansigmas[i,1,1]=mean(gens$gensigmas[(burnin+1):L,i][[1]][1,1]);
  meansigmas[i,1,2]=mean(gens$gensigmas[(burnin+1):L,i][[1]][1,2]);
  meansigmas[i,2,1]=mean(gens$gensigmas[(burnin+1):L,i][[1]][2,1]);
  meansigmas[i,2,2]=mean(gens$gensigmas[(burnin+1):L,i][[1]][2,2]);
}
# x[1,1,1]
#  mix1 <- normmix(ps=c(.3, .7), mus=list(c(0, 0), c(1, 1)),sigmas=list(.01*diag(2), .01*diag(2)))
#cat("\n data avg \n",mean(data[,1]),mean(data[,2]))
#windows()
for (i in 1:m)
{ 
  cat("\n posterior meanp\n",meanp[i])
  cat("\n trueps\n",truemix[[i]]$p)
  cat("\n posterior meanmus[",i,"]\n",meanmus[i,1],meanmus[i,2])
  cat("\n truemean[",i,"]\n",truemix[[i]]$mu)
  cat("\n posterior mean sigma[",i,"]\n",meansigmas[i,,])
  cat("\n true sigma[",i,"]\n",truemix[[i]]$sigma)
}
par(mfrow=c(1,2))
plot(gens$genmus[1,2,],type="l")
plot(gens$genmus[1,1,],type="l")
#}

}

GenNormalMixture<- function(n=50,m=1,  
                            xlims=c(-3,3),ylims=c(-3,3),r=20,trunc=FALSE)
{
  #  mix=GenNormalMixture(n,m,xlims=c(-10,10),ylims=c(-10,10),r=2,trunc)
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
  data<<-genNormMix_sppmix(n,mix1)
  opar <- par()      # make a copy of current settings
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
  par(opar)          
  return(mix1)
}

PlotNormalMixture<- function(mix1,data,
                             m=1,  
                             xlims=c(-3,3),
                             ylims=c(-3,3))
{
  #  PlotNormalMixture(mix1,data,m,xlims=c(-10,10),ylims=c(-10,10))
  n=nrow(data) 
  opar <- par()      # make a copy of current settings
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
  par(opar)          
}