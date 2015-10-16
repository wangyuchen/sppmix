//#include "sppmix.h"
#include <RcppArmadillo.h>
//#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <stdio.h>
#include <time.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat invmat2d(arma::mat const& A){
  arma::mat B=zeros(2,2);
  float det=A(0,0)*A(1,1)-A(0,1)*A(1,0);
  //  if(det<=0){Rcout << "\n"<<"A is not pd " << std::endl;}
  B(0,0)=A(0,0)/det;
  B(0,1)=-A(0,1)/det;
  B(1,0)=-A(1,0)/det;
  B(1,1)=A(1,1)/det;
  return B;
}

// [[Rcpp::export]]
arma::mat rnorm2(int n, arma::vec mu,
                 arma::mat sigma) {
  arma::mat Z = arma::randn(n, 2);
  double sig1=sqrt(sigma(0,0)),
    sig2=sqrt(sigma(1,1));
  double rho=sigma(0,1)/(sig1*sig2);
  //Rcout << " " << sigma << std::endl ;
  //gen bivariate normal, mu, sig
  arma::mat Gens(n, 2);
  colvec Z1 = Z.col(0),Z2 = Z.col(1);
  Gens.col(0) = sig1*Z1+mu(0);
  Gens.col(1) = sig2*(rho*Z1+
    sqrt(1-rho*rho)*Z2)+mu(1);
  //  if(det(sigma)<=0){Rcout << "\n"<<"sigma not pd "<<det(sigma)<<"\n"<<std::endl;}
  return Gens;
}

// [[Rcpp::export]]
arma::vec rDirichlet(arma::vec const& d){
  int k = d.size();
  arma::vec gens = zeros(k,1);
  for(int i=0;i<k;i++)
    gens(i)=rgamma(1,1/d(i))[0];
  return gens/sum(gens);
}

// [[Rcpp::export]]
arma::mat rWishart(int const& df, arma::mat const& A){
  arma::vec mu= zeros(2,1);
  arma::mat Gens=rnorm2(df, mu,A);
  return Gens.t()*Gens;
}

// [[Rcpp::export]]
int rBinom(int const& n,double const& p){
  arma::vec probs=zeros(n+1,1);
  double u=runif(1)[0],sum1=0;
  probs(0)=pow(1-p,n);
  sum1=probs(0);
  if (u<sum1)
    return 0;
  int gen1=0;
  for(int i=0;i<n;i++)
  {
    probs(i+1)=((n-i)/(i+1))*(p/(1-p))*probs(i);
    if (u>=sum1 && u<sum1+probs(i+1))
    {
      gen1=i+1;
      break;
    }
  }
  return gen1;
}

// [[Rcpp::export]]
double ApproxCompMass(int const& LL,arma::vec const& ticsx,
                      arma::vec const& ticsy,arma::mat const& areas,
                      arma::vec const& mu,
                      arma::mat const& sig,arma::mat const& siginv){
  double quad,approx=0,
    d1=1/sqrt(det(2*datum::pi*sig));
  int ii,jj;
  arma::vec midtics=zeros(2,1);
  for(jj=0;jj<LL-1;jj++)
  {
    for (ii=0;ii<LL-1;ii++)
    {
      midtics(0)=(ticsx(ii)+ticsx(ii+1))/2;
      midtics(1)=(ticsy(jj)+ticsy(jj+1))/2;
      quad=as_scalar((midtics-mu).t()*siginv*(midtics-mu));
      approx=approx+areas(ii,jj)*d1* exp(-.5*quad);
    }
  }
  return approx;
}

// [[Rcpp::export]]
double ApproxMHRatiomu(int const& LL,arma::vec const& ticsx,
                       arma::vec const& ticsy,arma::mat const& areas,
                       arma::vec const& curmu,arma::vec const& propmu,
                       arma::mat const& sig,arma::mat const& siginv){
  double quad,approxFmu=0,approxFpropmu=0,
    d1=1/sqrt(det(2*datum::pi*sig));
  int ii,jj;
  arma::vec midtics=zeros(2,1);
  for(jj=0;jj<LL-1;jj++)
  {
    for (ii=0;ii<LL-1;ii++)
    {
      midtics(0)=(ticsx(ii)+ticsx(ii+1))/2;
      midtics(1)=(ticsy(jj)+ticsy(jj+1))/2;
      quad=as_scalar((midtics-curmu).t()*siginv*(midtics-curmu));
      approxFmu=approxFmu+areas(ii,jj)*d1* exp(-.5*quad);
      quad=as_scalar((midtics-propmu).t()*siginv*(midtics-propmu));
      approxFpropmu=approxFpropmu+areas(ii,jj)*d1* exp(-.5*quad);
    }
  }
  return approxFmu/approxFpropmu;
}

// [[Rcpp::export]]
double ApproxMHRatiosig(int const& LL,arma::vec const& ticsx,
                        arma::vec const& ticsy,arma::mat const& areas,
                        arma::vec const& mu1,arma::mat const& propsigma,
                        arma::mat const& sig,arma::mat const& siginv){
  double quad,approxFsig=0,approxFpropsig=0,
    d1=1/sqrt(det(2*datum::pi*sig)),
    d1prop=1/sqrt(det(2*datum::pi*propsigma));
  int ii,jj;
  arma::mat invpropsigma=invmat2d(propsigma);
  arma::vec midtics=zeros(2,1);
  for(jj=0;jj<LL-1;jj++)
  {
    for (ii=0;ii<LL-1;ii++)
    {
      midtics(0)=(ticsx(ii)+ticsx(ii+1))/2;
      midtics(1)=(ticsy(jj)+ticsy(jj+1))/2;
      quad=as_scalar((midtics-mu1).t()*siginv*(midtics-mu1));
      approxFsig=approxFsig+areas(ii,jj)*d1* exp(-.5*quad);
      quad=as_scalar((midtics-mu1).t()*invpropsigma*(midtics-mu1));
      approxFpropsig=approxFpropsig+areas(ii,jj)*d1prop* exp(-.5*quad);
    }
  }
  return approxFsig/approxFpropsig;
}

// [[Rcpp::export]]
arma::vec rMultinomial(int const& n,arma::vec const& p){
  int j,k=p.size();
  arma::vec gen=zeros(k,1);
  /*  double S=1,l=n;
  for (j=0;j<k;j++)
  {
  if(p(j)<S)
  gen(j)=rBinom(l,p(j)/S);
  else
  gen(j)=1;
  if(gen(j)==1)
  break;
  else
  {
  l=l-gen(j);
  if(l==0)break;
  S=S-p(j);
  if(S<0)S=0;
  }
  }*/
  gen(0)=rBinom(n,p(0));
  //  Rcout << "\n"<<p(0)<<","<<gen(0) << std::endl ;
  int sum1=0;
  for(j=1;j<k-1;j++)
  {
    sum1=sum1+gen(j-1);
    gen(j)=rBinom(n-sum1,p(j)/sum(p.subvec(j,k-2)));
  }
  gen(k-1)=n-sum(gen);
  //      Rcout << " passed" << std::endl ;*/
  return gen;
}


//data, x-y limits, m=num of comps to fit
//L=num iter, burnin
//function [marginal,meanlamda,meanmus,meansigmas,
//meanp,mus,sigmas,ps,meanz,sumz,postmeanintensity,
//approx,lamdas,invsigs]=DAMCMC2dNormalMixture(
//data,xlimits,ylimits,m,L,burnin,trueintensity,
//trueps,truemus,truesigmas,truncate,skiplots)
//just get realizations no plotting here

//[[Rcpp::export]]
List DAMCMC2d(arma::mat const& data, arma::vec const& xlims,
              arma::vec const& ylims,
              int const& m,int const& L,
              int const& burnin,bool const& truncate)
{
  int n=data.n_rows;
  Rcout << "Dataset has " << n <<" points" << std::endl ;
  //  cube(n_rows, n_cols, n_slices)
  cube genz1=zeros(n,m,L),
    genmus=zeros(m,2,L);
  //    Rcout << " passed" << std::endl ;
  //  genmus.slice(0), first realization
  //all elements in the field are 2x2 matrices
  field<arma::mat> gensigmas(L,m),geninvsigmas(L,m);
  //  gensigmas(0,1), first realization for second comp
  arma::mat genz=zeros(n,m),
    genps=zeros(L,m),
    consts=zeros(L,m);

  int i,j,r,dat;
  ivec which=randi(m, distr_param(0,n-1));
  //process data for hyperparam values
  arma::vec mins=min(data,1);
  arma::vec maxs=max(data,1);
  double Rx=maxs(0)-mins(0),
    Ry=maxs(1)-mins(1);
  arma::vec ksi=zeros(2,1);
  ksi(0)=sum(data.col(0))/n;
  ksi(1)=sum(data.col(1))/n;
  //  Rcout << ksi<< std::endl ;
  arma::mat kappa,kappainv;
  kappa.eye(2,2);
  kappa(0,0)=1/(Rx*Rx);
  kappa(1,1)=1/(Ry*Ry);
  kappainv.eye(2,2);
  kappainv(0,0)=Rx*Rx;
  kappainv(1,1)=Ry*Ry;
  //hypers 3:a=3, 4:g=.3, 5:gam=1
  double a=3,g=1,gam=1;
  arma::mat hmat,hmatinv;
  hmat.eye(2,2);
  hmat(0,0)=100*g/(a*Rx*Rx);
  hmat(1,1)=100*g/(a*Ry*Ry);
  hmatinv.eye(2,2);
  hmatinv(0,0)=a*Rx*Rx/(100*g);
  hmatinv(1,1)=a*Ry*Ry/(100*g);
  //starting values
  for (i=0;i<m;i++)
  {
    //    genmus.slice(0);
    genmus(i,0,0)=data(which(i),0);
    genmus(i,1,0)=data(which(i),1);
    gensigmas(0,i)=kappainv;
    geninvsigmas(0,i)=invmat2d(gensigmas(0,i));
    genps(0,i)=1.0/m;
    consts(0,i)=1/sqrt(det(2*datum::pi*gensigmas(0,i)));
  }

  arma::mat sumsig,newdata,propmu=zeros(1,2),genmutemp=zeros(1,2),
    mutemp1=zeros(1,2),propsigma,
    zmultinomial=ones(n,m),sumxmu,ps1,
    ps2,beta,cov1;
  arma::vec qij=zeros(m,1),approx=zeros(m,1),
    ds,mu1=zeros(2,1),newmu=zeros(2,1);

  double ratio=1,approxcompj,sum1,quad;
  //setup grid for truncation
  int LL=101;
  arma::vec ticsx=zeros(LL),ticsy=zeros(LL);
  for (i=0;i<LL;i++)
  {
    ticsx(i)=mins(0)+i*(maxs(0)-mins(0))/(LL-1);
    ticsy(i)=mins(1)+i*(maxs(1)-mins(1))/(LL-1);
  }
  arma::mat areas=zeros(LL,LL);
  for(j=0;j<LL-1;j++)
    for(i=0;i<LL-1;i++)
      areas(i,j)=(ticsx(i+1)-ticsx(i))*(ticsy(j+1)-ticsy(j));

  //start MCMC
  Rcout << "Preliminaries done. Starting MCMC" << std::endl ;
  for(i=0;i<L-2;i++)
  {
    printf("\rWorking: %3.1f%% complete",100.0*i/(L-2));
    //  sample B matrix
    sumsig=zeros(2,2);
    for (j=0;j<m;j++)
      sumsig=sumsig+geninvsigmas(i,j);
    ps1=invmat2d(2*hmat+2*sumsig);
    //   if(det(ps1)<=0){Rcout << "\n"<<"ps1 not pd "<<det(ps1)<<"\n"<<std::endl;}
    beta=rWishart(2*g+2*m*a,ps1);
    //samples mus and sigmas
    ds=zeros(m,1);
    for(j=0;j<m;j++)
    {
      //samples mus
      sum1=sum(zmultinomial.col(j));
      //      Rcout <<"\n"<< i << " "<< j << " "<<sum1 << std::endl ;
      //     indi=find(zmultinomial.col(j)==1);
      if(sum1>0)
      {
        newdata=data.rows(find(zmultinomial.col(j)==1));
        //        Rcout << newdata<< std::endl ;
        //    if(newdata.n_rows!=sum1){Rcout << "\n"<<"not reading data properly "<<newdata.n_rows<< " "<<sum1 << std::endl ;}
        newmu(0)=sum(newdata.col(0))/sum1;
        newmu(1)=sum(newdata.col(1))/sum1;
        //        Rcout << " passeddfddd" << std::endl ;
      }
      else
      {
        newmu=zeros(2,1);
        //        Rcout << "\n"<<"Component with less than 2 points, exiting" << std::endl ;
        //        return List::create();
      }
      cov1=invmat2d(sum1*geninvsigmas(i,j)+kappa);
      //    if(det(cov1)<=0){Rcout << "\n"<<"Cov1 not pd, "<<cov1 << std::endl ;}
      //      Rcout << newmu << sum1<< std::endl ;
      mu1=cov1*(sum1*geninvsigmas(i,j)*newmu+kappa*ksi);
      //     Rcout << mu1<<cov1<< std::endl ;
      //proposed mu
      genmutemp=rnorm2(1,mu1,cov1);
      //      Rcout << " passed1" << std::endl ;
      if(truncate)
      {
        mu1(0)=genmus(j,0,i);
        mu1(1)=genmus(j,1,i);
        ratio=pow(ApproxMHRatiomu(LL,ticsx,ticsy,areas,
                                  mu1,trans(genmutemp),gensigmas(i,j),
                                  geninvsigmas(i,j)),sum1);
      }
      else
        ratio=1;
      if(runif(1)[0]<ratio)
      {
        genmus(j,0,i+1)=genmutemp(0);
        genmus(j,1,i+1)=genmutemp(1);
      }
      else
      {
        genmus(j,0,i+1)=genmus(j,0,i);
        genmus(j,1,i+1)=genmus(j,1,i);
      }
      //samples sigmas
      sumxmu=zeros(2,2);
      if (sum1>0)
        for(r=0;r<sum1;r++)
          sumxmu=sumxmu+(newdata.row(r)-genmutemp).t()*(newdata.row(r)-genmutemp);
      ps2=invmat2d(2*beta+sumxmu);
      //      if(det(ps2)<=0){Rcout << "\n"<<        "ps2 not pd "<<det(ps2)<<"\n"<<        det(beta)<< " "<<det(sumxmu)<<std::endl;}
      propsigma=rWishart(2*a+sum1,ps2);
      if(truncate)
      {
        mu1(0)=genmus(j,0,i+1);
        mu1(1)=genmus(j,1,i+1);
        ratio=pow(ApproxMHRatiosig(LL,ticsx,ticsy,areas,mu1,
                                   propsigma,gensigmas(i,j),
                                   geninvsigmas(i,j)),sum1);
      }
      else
        ratio=1;
      if(runif(1)[0]<ratio)
      {
        geninvsigmas(i+1,j)=propsigma;
        gensigmas(i+1,j)=invmat2d(geninvsigmas(i+1,j));
      }
      else
      {
        geninvsigmas(i+1,j)=geninvsigmas(i,j);
        gensigmas(i+1,j)=gensigmas(i,j);
      }
      ds(j)=gam+sum1;
      if (truncate)
      {
        mu1(0)=genmus(j,0,i+1);
        mu1(1)=genmus(j,1,i+1);
        approxcompj=ApproxCompMass(LL,ticsx,ticsy,
                                   areas,mu1,gensigmas(i+1,j),
                                   geninvsigmas(i+1,j));
      }
      else
        approxcompj=1;
      consts(i+1,j)=1/(approxcompj*
        sqrt(det(2*datum::pi*gensigmas(i+1,j))));
    }
    //sample component probs
    genps.row(i+1)=rDirichlet(ds).t();
    //sample indicators zij
    for(dat=0;dat<n;dat++)
    {
      for(j=0;j<m;j++)
      {
        genmutemp(0)=genmus(j,0,i+1);
        genmutemp(1)=genmus(j,1,i+1);
        //        Rcout << genmutemp<< std::endl ;
        //        Rcout << data.row(dat)<< std::endl ;
        mutemp1=data.row(dat)-genmutemp;
        quad=as_scalar(mutemp1*geninvsigmas(i+1,j)*trans(mutemp1));
        //        Rcout << mutemp1*geninvsigmas(i+1,j)*trans(mutemp1)<< std::endl ;
        qij(j)=genps(i+1,j)*consts(i+1,j)*exp(-.5*quad);
      }
      //      qij=qij/sum(qij);
      //      Rcout << qij<< std::endl ;
      //    Rcout << rMultinomial(1,qij)<< std::endl ;
      zmultinomial.row(dat)=reshape(rMultinomial(1,qij/sum(qij)),1,m);
    }
  }
  printf("\rDone                                                      ");

  return List::create(
    Named("genps") = genps,
    Named("genmus") = genmus,
    Named("gensigmas") = gensigmas,
    Named("genz") = genz);
}



/*** R
#major problem: need to allocate huge arrays
#this does not work, gives out of memory error
#Sys.setenv("PKG_CXXFLAGS"="-I../inst/include -DARMA_64BIT_WORD")
#Sys.setenv("PKG_LIBS"="$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)")
#Sys.setenv("CXX_STD"="CXX11")
# we're fine if we dont report all the z's

#library(Rcpp)
library(microbenchmark)
  library(MASS)
  mu=c(0,0)
  sig=rbind(c(5,0),c(0,2))
#rnorm_cpp(n,mu,sig)
  m <- diag(2);
#matrix(rnorm(4), 2, 2)
#getEigenValues(m)
#mvrnormArma(10,c(0,0),diag(2))
# microbenchmark(DAMCMC2d(n,mu,sig))
#  data
xlims=c(-3,3)
  ylims=c(-3,3)
  L=5000
m=1
burnin=1000
n=50
data=rnorm2(n,mu,sig)
# gens=DAMCMC2d(data,xlims,ylims,m,L,burnin,FALSE)
  meansmix(n,m,L,burnin,trunc=TRUE)
#  meansmix(n,m,L,burnin,trunc=FALSE)
  ***/
