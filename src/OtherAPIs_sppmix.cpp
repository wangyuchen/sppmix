//Written by Yuchen Wang, 2015
#include "sppmix.h"
//#include <Rcpp.h>
#include <mvtnormAPI.h>

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
double ApproxBivNormProb_sppmix(vec const& xlims,
                                vec const& ylims,vec const& mu,
                                mat const& sigma,int type)
{
  // NumericVector lls,
  //  NumericVector uls,NumericVector mu,
  //  NumericMatrix sigma
  NumericVector lls(2),uls(2);
  lls(0)=xlims(0);
  lls(1)=ylims(0);
  uls(0)=xlims(1);
  uls(1)=ylims(1);
  NumericVector muxy(2);
  muxy(0)=mu(0);
  muxy(1)=mu(1);
/*  NumericMatrix sig(2,2);
  sig(0,0)=sigma(0,0);
  sig(0,1)=sigma(0,1);
  sig(1,0)=sigma(1,0);
  sig(1,1)=sigma(1,1);*/
  int n = 2, nu = 0, maxpts = 2000, inform;
  int infin[2]={type,type};
  //INFIN INTEGER, array of integration limits flags:
  //if INFIN(I) < 0, Ith limits are (-infinity, infinity);
  //if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
  //if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
  //if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
  double abseps=1/1000, releps=1/1000, error, value;
  int rnd=1;
  double corr=sigma(0, 1) / sqrt(sigma(0, 0) * sigma(1, 1));
  /* mvtnorm_C_mvtdst is defined in mvtnorm/inst/include/mvtnormAPI.h */
  mvtnorm_C_mvtdst(&n, &nu, lls.begin(), uls.begin(), infin, &corr, muxy.begin(),
                   &maxpts, &abseps, &releps, &error, &value, &inform, &rnd);
  return value;
}
