#include "sppmix.h"

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

double densNormMixatx_sppmix(vec const& atx,
                             List const& mix)
{
  //  densNormMix_sppmix(c(0,0),truemix)
  int j,m=mix.size();
  List mth_comp;
  double c1,val=0,psj;
  vec muk,mu1=zeros(2);
  mat invsigk,sigmak;
  //  Rcout << m<< std::endl ;
  //  Rcout << atx<< std::endl ;
  for(j=0;j<m;j++)
  {
    mth_comp = mix[j];
    psj=as<double>(mth_comp["p"]);
    muk = as<vec>(mth_comp["mu"]);
    //    Rcout << j<< std::endl ;
    //    Rcout << muk<< std::endl ;
    sigmak = as<mat>(mth_comp["sigma"]);
    //    Rcout << sigmak<< std::endl ;
    mu1(0)=atx(0)-muk(0);
    mu1(1)=atx(1)-muk(1);
    //    Rcout << "passed"<< std::endl ;
    invsigk=invmat2d_sppmix(sigmak);
    c1=1.0/sqrt(det(2*datum::pi*sigmak));
    val=val+psj*c1*
      exp(-.5*as_scalar(trans(mu1)*invsigk*mu1));
  }
  return val;
}


mat dNormMix_sppmix(List const& mix, vec const& x,
                    vec const& y)
{
  int xnum=x.size(),ynum=y.size();
  //  int m = mix.size();
  //  List mth_comp;
  int i,j;
  mat z=zeros(xnum,ynum);
  vec atxy=zeros(2);
  for(i=0;i<xnum;i++)
    for(j=0;j<ynum;j++)
    {
      atxy(0)=x(i);
      atxy(1)=y(j);
      z(i,j)=densNormMixatx_sppmix(atxy,mix);
    }
    return z;
}

List GetStats_sppmix(vec const& gens,
                     double const& alpha)
{
  double mu=mean(gens);
  int L=gens.size();
  vec CS(2);
  vec sortedgens=sort(gens);
  //  Rcout << sortedgens<< std::endl ;
  CS(1)=sortedgens(floor((1-alpha/2)*L));
  CS(0)=sortedgens(floor(alpha*L/2));
  return List::create(
    Named("Min") = gens.min(),
    Named("Mean") = mu,
    Named("Max") = gens.max(),
    Named("CredibleSet") = CS,
    Named("CredibleSetConfidence") = 100*(1-alpha));
}

mat ApproxAvgPostIntensity(List const& genmix,
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
        intensityatxy=densNormMixatx_sppmix(xy,genmix[i]);
        AvgPostIntensity(x1,y1)=AvgPostIntensity(x1,y1)+intensityatxy/(L-burnin);
      }
      //      AvgPostIntensity(x1,y1)=AvgPostIntensity(x1,y1)/(L-burnin);
      countiter++;
    }
    return AvgPostIntensity;
  //,Named("PostIntensityofAvg") = PostIntensityAvg);
}


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

