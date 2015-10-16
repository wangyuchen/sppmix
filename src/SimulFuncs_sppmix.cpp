#include "sppmix.h"
//contains simulation functions only

// [[Rcpp::export]]
double rUnif_sppmix()
{
  // double u=randu<vec>(1)[0];
  //  return randu<vec>(1)[0];
  return ((double) rand() / (RAND_MAX+1.0));
}

// [[Rcpp::export]]
double rUnifab_sppmix(double const& a,
                      double const& b)
{
  // double u=randu<vec>(1)[0];
  return a+(b-a)*rUnif_sppmix();
}

// [[Rcpp::export]]
mat rnorm2_sppmix(int n,vec mu,mat sigma) {
 mat Z = randn(n, 2);
  double sig1=sqrt(sigma(0,0)),
    sig2=sqrt(sigma(1,1));
  double rho=sigma(0,1)/(sig1*sig2);
  //Rcout << " " << sigma << std::endl ;
  //gen bivariate normal, mu, sig
  mat Gens(n, 2);
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
  //      Rcout << " passed" << std::endl ;
  return gen;
}

// [[Rcpp::export]]
List rNormMix_sppmix(int const& lamda,
                       List const& mix)
{
  //mix[[k]] is the probability of a comp
  //mix[[k]] is a vec for mu, 2x1
  //mix[[k]] is a mat for sigma is 2x2
  int i,j,m  = mix.size(),comp=0;
  List mth_comp;
  //pick a component
  double u,sump,psj;
  int n=rpois(1,lamda)[0];
  //  Rcout << "mix" << std::endl ;
  //  Rcout << as<vec>(mth_comp["mu"]) << std::endl ;
  mat ret=zeros(n,2);
  vec comps(n);
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
    comps(i)=comp;
    //   Rcout << "comp="<<comp<< std::endl ;
    mth_comp = mix[comp];
    vec muk = as<vec>(mth_comp["mu"]);
    //  Rcout << muk<< std::endl ;
    mat sigmak = as<mat>(mth_comp["sigma"]);
    //  Rcout << sigmak<< std::endl ;
    ret.row(i)=rnorm2_sppmix(1,muk,sigmak);
  }
  return List::create(
    Named("data") = ret,
    Named("comp") = comps);
}

// [[Rcpp::export]]
vec rPerm_sppmix(int const& n)
{
  vec perm(n);
  for (int i=0; i<n; i++) perm(i)=i+1;
  std::random_shuffle(perm.begin(),perm.end());
  return perm;
}
