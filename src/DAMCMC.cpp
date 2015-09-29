//#include "sppmix.h"
#include <RcppArmadillo.h>
//#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <stdio.h>
#include <time.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
mat invmat2d_sppmix(mat const& A){
  mat B=zeros(2,2);
  double det=A(0,0)*A(1,1)-A(0,1)*A(1,0);
  //  if(det<=0){Rcout << "\n"<<"A is not pd " << std::endl;}
  B(0,0)=A(1,1)/det;
  B(0,1)=-A(0,1)/det;
  B(1,0)=-A(1,0)/det;
  B(1,1)=A(0,0)/det;
  return B;
}

// [[Rcpp::export]]
double rUnifab_sppmix(double const& a,
                    double const& b)
{
 // double u=randu<vec>(1)[0];
  return a+(b-a)*randu<vec>(1)[0];
}

// [[Rcpp::export]]
double rUnif_sppmix()
{
  // double u=randu<vec>(1)[0];
  return randu<vec>(1)[0];
}

// [[Rcpp::export]]
arma::mat rnorm2_sppmix(int n, arma::vec mu,
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
mat rWishart_sppmix(int const& df, mat const& A){
  mat Gens=rnorm2_sppmix(df, zeros(2),A);
  return Gens.t()*Gens;
}

// [[Rcpp::export]]
int rDiscrete_sppmix(int const& start,
                     int const& end,
                     vec const& probs){
  double u=rUnif_sppmix(), cdf=0;
  int gen1=start,m=probs.size();
  for(int i=0;i<m;i++)
  {
    if (u>cdf && u<=cdf+probs(i))
    {
      gen1=start+i;
      break;
    }
    cdf=cdf+probs(i);
  }
  return gen1;
}

// [[Rcpp::export]]
int rBinom_sppmix(int const& n,
                  double const& p){
  vec probs(n+1);
  probs(0)=pow(1-p,n);
  for(int i=1;i<n;i++)
    probs(i)=((n-i+1.0)/i)*(p/(1-p))*probs(i-1);
  probs(n)=pow(p,n);
//  Rcout << probs<< std::endl ;
  return rDiscrete_sppmix(0,n,probs);
}

// [[Rcpp::export]]
double rGamma_sppmix(double const& a,
                  double const& b){
  double u=rUnif_sppmix();
  double y=- a * b * log(rUnif_sppmix());
  double u0=pow(y,a - 1)*
    exp(-y*(a - 1)/(a * b))
            *exp(a - 1)/pow(a*b,a-1);
    while(u >= u0) {
      y=- a * b * log(rUnif_sppmix());
      u0=pow(y,a - 1)*
        exp(-y*(a - 1)/(a * b))
        *exp(a - 1)/pow(a*b,a-1);
      u=rUnif_sppmix();
    }
  return y;
}

// [[Rcpp::export]]
double rExp_sppmix(double const& a)
{
  return -log(rUnif_sppmix())/a;
}

// [[Rcpp::export]]
vec rDirichlet_sppmix(vec const& d){
  int k = d.size();
  vec gens(k);
  for(int i=0;i<k;i++)
    gens(i)=//rgamma(1.0,1.0/d(i))[0];
      rExp_sppmix(1.0/d(i));
  return gens/sum(gens);
}

// [[Rcpp::export]]
mat genNormMix_sppmix(int const& n,List const& mix)
  {
  //mix[[k]] is the probability of a comp
  //mix[[k]] is a vec for mu, 2x1
  //mix[[k]] is a mat for sigma is 2x2
  int i,j,m  = mix.size(),comp=0;
  List mth_comp;
//pick a component
  double u,sump,psj;
//  Rcout << "mix" << std::endl ;
//  Rcout << as<vec>(mth_comp["mu"]) << std::endl ;
  mat ret=zeros(n,2);
  for(i=0;i<n;i++)
  {
    u=rUnif_sppmix();
    sump=0;
    for(j=0;j<m;j++)
    {
      mth_comp = mix[j];
      psj=as<double>(mth_comp["p"]);
  //    Rcout << "psj="<<psj<< std::endl ;
      if (sump<u && u<=sump+psj)
      {
        comp=j;
        break;
      }
      sump=sump+psj;
    }
   //   Rcout << "comp="<<comp<< std::endl ;
    mth_comp = mix[comp];
    vec muk = as<vec>(mth_comp["mu"]);
  //  Rcout << muk<< std::endl ;
    mat sigmak = as<mat>(mth_comp["sigma"]);
  //  Rcout << sigmak<< std::endl ;
    ret.row(i)=rnorm2_sppmix(1,muk,sigmak);
  }
  return ret;
  //   as<vec>(mth_comp["mu"]),
   //  as<mat>(mth_comp["sigma"]));
}

// [[Rcpp::export]]
double ApproxCompMass_sppmix(int const& LL,vec const& ticsx,
                         vec const& ticsy,mat const& areas,
                         vec const& mu,
                         mat const& sig,mat const& siginv){
    double quad,approx=0,
      d1=1/sqrt(det(2*datum::pi*sig));
    int ii,jj;
    vec midtics=zeros(2,1);
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
double ApproxMHRatiomu_sppmix(int const& LL,vec const& ticsx,
                vec const& ticsy,mat const& areas,
                vec const& curmu,vec const& propmu,
                mat const& sig,mat const& siginv){
  double quad,approxFmu=0,approxFpropmu=0,
    d1=1/sqrt(det(2*datum::pi*sig));
  int ii,jj;
  vec midtics=zeros(2,1);
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
double ApproxMHRatiosig_sppmix(int const& LL,vec const& ticsx,
                         vec const& ticsy,mat const& areas,
                         vec const& mu1,mat const& propsigma,
                         mat const& sig,mat const& siginv){
    double quad,approxFsig=0,approxFpropsig=0,
      d1=1/sqrt(det(2*datum::pi*sig)),
      d1prop=1/sqrt(det(2*datum::pi*propsigma));
    int ii,jj;
    mat invpropsigma=invmat2d_sppmix(propsigma);
    vec midtics=zeros(2,1);
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
ivec rMultinomial_sppmix(int const& n,vec const& ps){
  int j,i,k=ps.size();
  ivec gen = zeros<ivec>(k) ;
//  if(sum(p)<1)
//    Rcout << "\n"<<p<<","<<sum(p) << std::endl ;
  gen(0)=rBinom_sppmix(n,ps(0));
  int sum1=0;
  double sump;
  for(j=1;j<k-1;j++)
  {
    sum1=sum1+gen(j-1);
    if (sum1==n)break;
    sump=0;
    for(i=j;i<k;i++)sump=sump+ps(i);
    gen(j)=rBinom_sppmix(n-sum1,ps(j)/sump);
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
List DAMCMC2d_sppmix(mat const& data, 
              vec const& xlims,
              vec const& ylims,
              int const& m,int const& L,
              int const& burnin,bool const& truncate)
{
  int n=data.n_rows;
  Rcout << "Dataset has " << n <<" points" << std::endl ;
//  cube(n_rows, n_cols, n_slices)
  cube //genz1=zeros(n,m,L),
       genmus=zeros(m,2,L);
//    Rcout << " passed" << std::endl ;
  //  genmus.slice(0), first realization
  //all elements in the field are 2x2 matrices
  field<mat> gensigmas(L,m),geninvsigmas(L,m);
//  gensigmas(0,1), first realization for second comp
  mat genz=zeros(n,m),
      genps=zeros(L,m),
      consts=zeros(L,m);

  int i,j,r,dat;
  ivec which=randi(m, distr_param(0,n-1));
//process data for hyperparam values
  vec mins=zeros(2),maxs=zeros(2);
  if(xlims.size()==1)//passed 0, use data
  {
    mins(0)=min(data.col(0));
    mins(1)=min(data.col(1));
    maxs(0)=max(data.col(0));
    maxs(1)=max(data.col(1));
  }else
  {
    mins(0)=xlims(0);
    mins(1)=ylims(0);
    maxs(0)=xlims(1);
    maxs(1)=ylims(1);
  }
//  Rcout << "mins="<<mins<< std::endl ;
  double Rx=max(data.col(0))-min(data.col(0)),
         Ry=max(data.col(1))-min(data.col(1));
  vec ksi=zeros(2,1);
  ksi(0)=sum(data.col(0))/n;
  ksi(1)=sum(data.col(1))/n;
//  Rcout << ksi<< std::endl ;
  mat kappa,kappainv;
  kappa.eye(2,2);
  kappa(0,0)=1/(Rx*Rx);
  kappa(1,1)=1/(Ry*Ry);
  kappainv.eye(2,2);
  kappainv(0,0)=Rx*Rx;
  kappainv(1,1)=Ry*Ry;
 //hypers 3:a=3, 4:g=.3, 5:gam=1
  double a=3,g=1,gam=1;
  imat prevz,zmultinomial(n,m);
  mat hmat,hmatinv;
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
    geninvsigmas(0,i)=kappa;//invmat2d_sppmix(gensigmas(0,i));
    genps(0,i)=1.0/m;
    consts(0,i)=1.0/sqrt(det(2*datum::pi*gensigmas(0,i)));
  }
  for (dat=0;dat<n;dat++)
    zmultinomial.row(dat)=reshape(rMultinomial_sppmix(1,rDirichlet_sppmix(ones(m))),1,m);
//  Rcout << sum(zmultinomial.col(0))<< std::endl ;
//  Rcout << sum(zmultinomial.col(1))<< std::endl ;
  //Rcout << zmultinomial<< std::endl ;
//    return List::create();
  mat sumsig,newdata,propmu=zeros(1,2),
    genmutemp=zeros(1,2),
    mutemp1=zeros(1,2),propsigma,
    sumxmu,ps1,ps2,beta,cov1;
  vec qij(m),approx=zeros(m),
    ds,qs=zeros(m),mu1=zeros(2),
    newmu=zeros(2);
  
  double MHjump=0,ratio=1,sumd,approxcompj;//,quad;
//setup grid for truncation
  int LL=21,sum1;
  vec ticsx=zeros(LL),ticsy=zeros(LL);
  for (i=0;i<LL;i++)
  {
    ticsx(i)=mins(0)+i*(maxs(0)-mins(0))/(LL-1);
    ticsy(i)=mins(1)+i*(maxs(1)-mins(1))/(LL-1);
  }
  mat areas=zeros(LL,LL);
  for(j=0;j<LL-1;j++)
    for(i=0;i<LL-1;i++)
      areas(i,j)=(ticsx(i+1)-ticsx(i))*(ticsy(j+1)-ticsy(j));
 
//start MCMC
  Rcout << "Preliminaries done. Starting MCMC" << std::endl ;
  for(i=0;i<L-1;i++)
  {
    prevz=zmultinomial;
    printf("\rWorking: %3.1f%% complete (iteration %d/%d)",100.0*i/(L-2),i,L);
//  sample B matrix
    sumsig=zeros(2,2);
    for (j=0;j<m;j++)
      sumsig=sumsig+geninvsigmas(i,j);
    ps1=invmat2d_sppmix(2*hmat+2*sumsig);
//   if(det(ps1)<=0){Rcout << "\n"<<"ps1 not pd "<<det(ps1)<<"\n"<<std::endl;}
    beta=rWishart_sppmix(2*g+2*m*a,ps1);
//samples mus and sigmas
    ds=zeros(m,1);
    for(j=0;j<m;j++)
    {
//samples mus
      sum1=sum(zmultinomial.col(j));
      //     indi=find(zmultinomial.col(j)==1);
      if(sum1>0)
      {
        newdata=data.rows(find(zmultinomial.col(j)==1));
//        Rcout << newdata<< std::endl ;
//    if(newdata.n_rows!=sum1){Rcout << "\n"<<"not reading data properly "<<newdata.n_rows<< " "<<sum1 << std::endl ;return List::create();}
        newmu(0)=sum(newdata.col(0))/sum1;
        newmu(1)=sum(newdata.col(1))/sum1;
//        Rcout << newmu << std::endl ;
      }
      else
      {
 //       break;
// Rcout <<"\n"<< i << " "<< j << " "<<sum1 << std::endl ;
 //         return List::create();
        newmu=zeros(2,1);
//        Rcout << "\n"<<"Component with less than 2 points, exiting" << std::endl ;
//        return List::create();
      }
      cov1=invmat2d_sppmix(sum1*geninvsigmas(i,j)+kappa);
//    if(det(cov1)<=0){Rcout << "\n"<<"Cov1 not pd, "<<cov1 << std::endl ;}
    //      Rcout << newmu << sum1<< std::endl ;
      mu1=cov1*(sum1*geninvsigmas(i,j)*newmu+kappa*ksi);
//      Rcout << mu1<<cov1<< std::endl ;
 //proposed mu
      genmutemp=rnorm2_sppmix(1,mu1,cov1);
//      Rcout << "genmu" <<genmutemp<< std::endl ;
      if(truncate)
      {
        mu1(0)=genmus(j,0,i);
        mu1(1)=genmus(j,1,i);
        ratio=pow(ApproxMHRatiomu_sppmix(LL,ticsx,ticsy,areas,
          mu1,trans(genmutemp),gensigmas(i,j),
          geninvsigmas(i,j)),sum1);
      }
      else 
        ratio=1;
      if(rUnif_sppmix()<ratio)
      {
        genmus(j,0,i+1)=genmutemp(0);
        genmus(j,1,i+1)=genmutemp(1);
      }
      else
      {
        genmutemp(0)=genmus(j,0,i);
        genmutemp(1)=genmus(j,1,i);
        genmus(j,0,i+1)=genmus(j,0,i);
        genmus(j,1,i+1)=genmus(j,1,i);
      }
//samples sigmas
      sumxmu=zeros(2,2);
      if (sum1>0)
        for(r=0;r<sum1;r++)
          sumxmu=sumxmu+(newdata.row(r)-genmutemp).t()*(newdata.row(r)-genmutemp);
//      Rcout << sumxmu<<"\n"<<std::endl;
      ps2=invmat2d_sppmix(2*beta+sumxmu);
      cov1=rWishart_sppmix(2*a+sum1,ps2);
//      if(det(ps2)<=0){Rcout << "\n"<<        "ps2 not pd "<<det(ps2)<<"\n"<<        det(beta)<< " "<<det(sumxmu)<<std::endl;}
      propsigma=invmat2d_sppmix(cov1);
      if(truncate)
      {
        mu1(0)=genmus(j,0,i+1);
        mu1(1)=genmus(j,1,i+1);
        ratio=pow(ApproxMHRatiosig_sppmix(LL,ticsx,ticsy,areas,mu1,
                                   propsigma,gensigmas(i,j),
                                  geninvsigmas(i,j)),sum1);
      }
      else 
        ratio=1;
      if(rUnif_sppmix()<ratio)
      {
        geninvsigmas(i+1,j)=cov1;//invmat2d_sppmix(propsigma);
        gensigmas(i+1,j)=propsigma;
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
        approxcompj=ApproxCompMass_sppmix(LL,ticsx,ticsy,
                            areas,mu1,gensigmas(i+1,j),
                            geninvsigmas(i+1,j));
      }
      else
        approxcompj=1;
      consts(i+1,j)=1.0/(approxcompj*
          sqrt(det(2*datum::pi*gensigmas(i+1,j))));
    }
//sample component probs
    genps.row(i+1)=rDirichlet_sppmix(ds).t();
//    Rcout << genps.row(i+1)<< std::endl ;
//    Rcout << ds<< std::endl ;
    //   Rcout << "1 "<<sum(zmultinomial.col(0))<< std::endl ;
 //   Rcout << "2 "<<sum(zmultinomial.col(1))<< std::endl ;
//    Rcout << consts.row(i+1)<< std::endl ;
 
    //sample indicators zij
    for(dat=0;dat<n;dat++)
    {
      sumd=0;
      for(j=0;j<m;j++)
      {
//        genmutemp(0)=genmus(j,0,i+1);
 //       genmutemp(1)=genmus(j,1,i+1);
//        Rcout << genmutemp<< std::endl ;
//        Rcout << data.row(dat)<< std::endl ;
//        mutemp1=data.row(dat)-genmutemp;
        mutemp1(0)=data(dat,0)-genmus(j,0,i+1);
        mutemp1(1)=data(dat,1)-genmus(j,1,i+1);
        //        quad=as_scalar(mutemp1*geninvsigmas(i+1,j)*trans(mutemp1));
//        Rcout << quad<< std::endl ;
//        Rcout << mutemp1<< std::endl ;
        qij(j)=genps(i+1,j)*consts(i+1,j)*
          exp(-.5*as_scalar(mutemp1*geninvsigmas(i+1,j)*trans(mutemp1)));
        sumd=sumd+qij(j);
//        Rcout << qij(dat,j)<< std::endl ;
      }
      qs=qij/sumd;
//if(sum(qs)<1)
//Rcout << sum(qs)<< std::endl ;
//      if (sum(qij)==1)
      zmultinomial.row(dat)=reshape(rMultinomial_sppmix(1,qs),1,m);
//      zmultinomial.row(dat)=reshape(rMultinomial_sppmix(1,qij/sum(qij)),1,m);
//      Rcout << sum(zmultinomial.row(dat))<< std::endl ;
    }
    ratio=1;
    for(j=0;j<m;j++)
    {
      //component with less that 2 points is not accepted
      if (sum(zmultinomial.col(j))<1)
      {
//        Rcout <<i<<zmultinomial<< std::endl ;
 //       Rcout <<"\nComponent "<<j+1<<" has less than 2 pts"<< std::endl ;
 //       Rcout << "1 "<<sum(zmultinomial.col(0))<< std::endl ;
 //       Rcout << "2 "<<sum(zmultinomial.col(1))<< std::endl ;
 //       return List::create();
        ratio=0;
        break;
      }
    }
    if(rUnif_sppmix()<ratio)
    {
      MHjump=MHjump+1;
      genz=genz+zmultinomial;
    }
    else
    {
      genps.row(i+1)=genps.row(i);
      for(j=0;j<m;j++)
      {
        geninvsigmas(i+1,j)=geninvsigmas(i,j);
        gensigmas(i+1,j)=gensigmas(i,j);
        genmus(j,0,i+1)=genmus(j,0,i);
        genmus(j,1,i+1)=genmus(j,1,i);
      }
      genz=genz+prevz;
      zmultinomial=prevz;
    }
  }
  printf("\rDone                                                      ");
  printf("\rMH acceptance %3.1f%%",MHjump/L);
  genz=genz/L;
  
  return List::create(
    Named("genps") = genps,
    Named("genmus") = genmus,
    Named("gensigmas") = gensigmas,
    Named("genz") = genz,
    Named("ticsnum") = LL,
    Named("ticsx") = ticsx,
    Named("ticsy") = ticsy,
    Named("ticsareas") = areas);
}



/*** R
#major problem: need to allocate huge arrays
#this does not work, gives out of memory error
#Sys.setenv("PKG_CXXFLAGS"="-I../inst/include -DARMA_64BIT_WORD")
#Sys.setenv("PKG_LIBS"="$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)")
#Sys.setenv("CXX_STD"="CXX11")
# we're fine if we dont report all the z's
#stl


#RUN HELPERFUNCS.R FIRST

#library(Rcpp)
library(microbenchmark)
library(MASS)
  mu=c(0,0)
  sig=rbind(c(1,0),c(0,1))
#rnorm_cpp(n,mu,sig)
  m <- diag(2)
  #matrix(rnorm(4), 2, 2)
#getEigenValues(m)
#mvrnormArma(10,c(0,0),diag(2))
 # microbenchmark(DAMCMC2d(n,mu,sig))
#  data
  xlims<-c(0,10)
  ylims<-c(0,10)
  L=5000
  m=2
  burnin=1000
  n=100
 # data=rnorm2_sppmix(n,mu,sig)
  graphics.off()
  truemix=GenNormalMixture(n=100,m=2,xlims=c(0,10),ylims=c(0,10),r=2,trunc)
#  PlotNormalMixture(mix1,data,m,xlims=c(0,10),ylims=c(0,10))
  # gens=DAMCMC2d(data,xlims,ylims,m,L,burnin,FALSE)
  windows()
#  meansmix(truemix,n,m,xlims,ylims,L,burnin,trunc=TRUE)
  meansmix(truemix,n,m,xlims,ylims,L,burnin,trunc=FALSE)
  */
