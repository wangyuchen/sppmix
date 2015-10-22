#include "sppmix.h"
//approximation functions only
//used by cpp and R functions
//visible in the package only

//[[Rcpp::export]]
mat ApproxAvgPostIntensity(List const& genmix,
                           vec const& lamdas,
                           int const& LL,
                           int const& burnin,
                           vec const& ticsx,
                           vec const& ticsy)
{
  //compute the average of the posterior surfaces
  //needs to be multiplied by the lamdas
  int countiter=0,i,L  = genmix.size();
  mat AvgPostIntensity = zeros(LL,LL);
  //,PostIntensityAvg = zeros(LL,LL);
  vec xy(2);
  //  Rcout << muk<< std::endl ;
  //  List tmix=genmix[0];
  //  int m=tmix.size();
  //  Rcout <<m<< std::endl ;
  //  List mixcomp(m);//list containing mixture ps,mus,sigs
  //  List mixcomp;//list containing mixture ps,mus,sigs
  double intensityatxy;
  //  countiter=0;
  for(int x1=0;x1<LL;x1++)
    for(int y1=0;y1<LL;y1++)
    {
      printf("\rComputing intensity surfaces: %3.1f%% complete",100.0*countiter/(LL*LL));
      xy(0)=ticsx(x1);
      xy(1)=ticsy(y1);
      //for each realization, compute the intensity
      //surface at the point xy
      for(i=burnin;i<L;i++)
      {
        //        mixcomp=genmix[i];
        //        intensityatxy=densNormMixatx_sppmix(xy,mixcomp);
        intensityatxy=lamdas(i)*densNormMixatx_sppmix(xy,genmix[i]);
        AvgPostIntensity(x1,y1)=AvgPostIntensity(x1,y1)+intensityatxy/(L-burnin);
      }
      countiter++;
    }
    printf("\rDone                                                      \n");
  return AvgPostIntensity;
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
