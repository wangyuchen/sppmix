#' @export
meansmix<- function(n=50,m=1,L=10000,burnin=1000,trunc=FALSE)
{
  #  L=5000
  #  m=1
  #  burnin=1000
  #  n=50
  mu=c(0,0)
  sig=rbind(c(5,0),c(0,2))
  #rnorm_cpp(n,mu,sig)
  # m <- diag(2);
  #matrix(rnorm(4), 2, 2)
  #getEigenValues(m)
  #mvrnormArma(10,c(0,0),diag(2))
  # microbenchmark(DAMCMC2d(n,mu,sig))
  data=rnorm2(n,mu,sig)
  #  data
  xlims=c(-3,3)
  ylims=c(-3,3)
  gens=DAMCMC2d(data,xlims,ylims,m,L,burnin,trunc)

  meanmus=matrix(0,m,2)
  meansigmas=matrix(0,m,2,2)
  # meaninvsigmas=zeros(m,2,2);
  meanp=matrix(0,m,1)
  for (i in 1:m)
  {
    meanp[i]=mean(gens$genps[(burnin+1):L,i]);
    meanmus[i,1]=mean(gens$genmus[i,1,(burnin+1):L]);
    meanmus[i,2]=mean(gens$genmus[i,2,(burnin+1):L]);
    # meansigmas(i,1,1)=mean(gens$gensigmas[burnin+1:L,i]((burnin+1):L,i,1,1));
    #  meansigmas(i,1,2)=mean(gens$gensigmas((burnin+1):L,i,1,2));
    # meansigmas(i,2,1)=mean(gens$gensigmas((burnin+1):L,i,2,1));
    #  meansigmas(i,2,2)=mean(gens$gensigmas((burnin+1):L,i,2,2));
    # meansigma=reshape(meansigmas(i,:,:),2,2)
    #   meaninvsigmas(i,:,:)=meansigma^-1;
  }
  # x[1,1,1]
  #  mix1 <- normmix(ps=c(.3, .7), mus=list(c(0, 0), c(1, 1)),sigmas=list(.01*diag(2), .01*diag(2)))
  cat("\n meanp\n",meanp)
  cat("\n meanmus\n",meanmus[1,1],meanmus[1,2])
  cat("\n original means \n",mean(data[,1]),mean(data[,2]))
  par(mfrow=c(1,2))
  plot(gens$genmus[1,2,])
  plot(gens$genmus[1,1,])

}
