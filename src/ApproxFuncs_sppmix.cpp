#include "sppmix.h"
//Written by Sakis Micheas, 2015
//approximation functions only
//used by cpp and R functions
//visible in the package only

//' @export
//[[Rcpp::export]]
mat ApproxAvgPostIntensity(List const& genmix,
                           vec const& lamdas,
                           int const& LL,
                           int const& burnin,
                           vec const& xlims,
                           vec const& ylims)
{
  //compute the average of the posterior surfaces
  //needs to be multiplied by the lamdas
  int countiter=0,i,L  = genmix.size();
  mat AvgPostIntensity = zeros(LL,LL);
  //,PostIntensityAvg = zeros(LL,LL);
  vec xy(2);
  vec ticsx=zeros(LL),ticsy=zeros(LL);
  for (i=0;i<LL;i++)
  {
    ticsx(i)=xlims(0)+i*(xlims(1)-xlims(0))/(LL-1);
    ticsy(i)=ylims(0)+i*(ylims(1)-ylims(0))/(LL-1);
  }

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

//' @export
// [[Rcpp::export]]
double ApproxCompMass_sppmix(int const& LL,
  vec const& xlims,vec const& ylims,vec const& mu,
  mat const& sig,mat const& siginv)
{
//approximates a single multivariate normal mass
  double quad,approx=0,
    d1=1/sqrt(det(2*datum::pi*sig));
  int i,j;
  vec midtics=zeros(2,1);
  vec ticsx=zeros(LL),ticsy=zeros(LL);
  for (i=0;i<LL;i++)
  {
    ticsx(i)=xlims(0)+i*(xlims(1)-xlims(0))/(LL-1);
    ticsy(i)=ylims(0)+i*(ylims(1)-ylims(0))/(LL-1);
  }

  for(j=0;j<LL-1;j++)
  {
    for (i=0;i<LL-1;i++)
    {
      double area=(ticsx(i+1)-ticsx(i))*(ticsy(j+1)-ticsy(j));
      midtics(0)=(ticsx(i)+ticsx(i+1))/2;
      midtics(1)=(ticsy(j)+ticsy(j+1))/2;
      quad=as_scalar((midtics-mu).t()*siginv*(midtics-mu));
      approx=approx+area*d1* exp(-.5*quad);
    }
  }
  return approx;
}

//' @export
// [[Rcpp::export]]
double ApproxMHRatiomu_sppmix(int const& LL,
  vec const& xlims,vec const& ylims,
  vec const& curmu,vec const& propmu,
  mat const& sig,mat const& siginv)
{
  double quad,approxFmu=0,approxFpropmu=0,
    d1=1/sqrt(det(2*datum::pi*sig));
  int i,j;
  vec ticsx=zeros(LL),ticsy=zeros(LL);
  for (i=0;i<LL;i++)
  {
    ticsx(i)=xlims(0)+i*(xlims(1)-xlims(0))/(LL-1);
    ticsy(i)=ylims(0)+i*(ylims(1)-ylims(0))/(LL-1);
  }
  vec midtics=zeros(2,1);
  for(j=0;j<LL-1;j++)
  {
    for (i=0;i<LL-1;i++)
    {
      double area=(ticsx(i+1)-ticsx(i))*(ticsy(j+1)-ticsy(j));
      midtics(0)=(ticsx(i)+ticsx(i+1))/2;
      midtics(1)=(ticsy(j)+ticsy(j+1))/2;
      quad=as_scalar((midtics-curmu).t()*siginv*(midtics-curmu));
      approxFmu=approxFmu+area*d1* exp(-.5*quad);
      quad=as_scalar((midtics-propmu).t()*siginv*(midtics-propmu));
      approxFpropmu=approxFpropmu+area*d1* exp(-.5*quad);
    }
  }
  return approxFmu/approxFpropmu;
}

//' @export
// [[Rcpp::export]]
double ApproxMHRatiosig_sppmix(int const& LL,
  vec const& xlims,vec const& ylims,vec const& mu1,
  mat const& propsigma,mat const& sig,
  mat const& siginv)
{
  double quad,approxFsig=0,approxFpropsig=0,
    d1=1/sqrt(det(2*datum::pi*sig)),
    d1prop=1/sqrt(det(2*datum::pi*propsigma));
  mat invpropsigma=invmat2d_sppmix(propsigma);
  int i,j;
  vec ticsx=zeros(LL),ticsy=zeros(LL);
  for (i=0;i<LL;i++)
  {
    ticsx(i)=xlims(0)+i*(xlims(1)-xlims(0))/(LL-1);
    ticsy(i)=ylims(0)+i*(ylims(1)-ylims(0))/(LL-1);
  }
  vec midtics=zeros(2,1);
  for(j=0;j<LL-1;j++)
  {
    for (i=0;i<LL-1;i++)
    {
      double area=(ticsx(i+1)-ticsx(i))*(ticsy(j+1)-ticsy(j));
      midtics(0)=(ticsx(i)+ticsx(i+1))/2;
      midtics(1)=(ticsy(j)+ticsy(j+1))/2;
      quad=as_scalar((midtics-mu1).t()*siginv*(midtics-mu1));
      approxFsig=approxFsig+area*d1* exp(-.5*quad);
      quad=as_scalar((midtics-mu1).t()*invpropsigma*(midtics-mu1));
      approxFpropsig=approxFpropsig+area*d1prop* exp(-.5*quad);
    }
  }
  return approxFsig/approxFpropsig;
}


//' @export
//[[Rcpp::export]]
mat ApproxBayesianModelAvgIntensity_sppmix(
    List const& genBDmix,vec const& lamdas,
    vec const& numcomp,vec const& distr_numcomp,
    int const& mincomp,int const& maxcomp,
    int const& LL,vec const& xlims,vec const& ylims)
{
  //apply burnin before calling this function
  //compute the average of the posterior surfaces
  //needs to be multiplied by the lamdas
  //and weighted by the comp number relative frequency
  //mincomp, maxcomp are 1,2,3,...,maxnumcomp, integers
  int countiter=0,i;
  mat AvgPostIntensity = zeros(LL,LL);
  vec ticsx=zeros(LL),ticsy=zeros(LL);
  for (i=0;i<LL;i++)
  {
    ticsx(i)=xlims(0)+i*(xlims(1)-xlims(0))/(LL-1);
    ticsy(i)=ylims(0)+i*(ylims(1)-ylims(0))/(LL-1);
  }
  vec xy(2);
  double weight,intensityatxy;
  printf("\nComputing Bayesian model average surface\n");
  for(int kval=mincomp-1;kval<maxcomp;kval++)
  {
    mat PostIntensity = zeros(LL,LL);
    uvec indi=find(numcomp==kval);
    if(sum(indi)==0)//no realizations for this k
      continue;
    int nn1=indi.size();
    vec newlamdas=lamdas(indi);
    List newgenmix(nn1);
    for(i=0;i<nn1;i++)
      newgenmix[i]=genBDmix[indi(i)];
    for(int x1=0;x1<LL;x1++)
    {
      for(int y1=0;y1<LL;y1++)
      {
        printf("\rNumber of components=%2d, surface count=%5d, %3.1f%% complete",kval,nn1,100.0*countiter/(LL*LL*(maxcomp-mincomp+1)));
//        printf("\rComputing Bayesian model average surface: %3.1f%% complete",100.0*countiter/(LL*LL*(maxcomp-mincomp+1)));
        xy(0)=ticsx(x1);
        xy(1)=ticsy(y1);
        //for each realization, compute the intensity
        //surface at the point xy
        for(i=0;i<nn1;i++)
        {
          intensityatxy=newlamdas(i)*densNormMixatx_sppmix(xy,newgenmix[i]);
          PostIntensity(x1,y1)=PostIntensity(x1,y1)+intensityatxy/nn1;
        }
        countiter++;
      }
    }
    weight=distr_numcomp(kval);
    AvgPostIntensity=AvgPostIntensity+weight*PostIntensity;
  }
  printf("\rDone                                                                     \n");
  return AvgPostIntensity;
}
