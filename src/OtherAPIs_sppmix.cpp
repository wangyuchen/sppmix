//Written by Yuchen Wang, 2015
#include <Rcpp.h>
#include <mvtnormAPI.h>

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
double pmvnorm_mvtnorm(NumericVector lls, NumericVector uls,
                       NumericVector mu, NumericMatrix sigma)
{

  int n = 2, nu = 0, maxpts = 2000, inform;
  int infin [2] = {2, 2};
  double abseps = 1/1000, releps = 1/1000, error, value;
  int rnd = 1;

  double corr = sigma(0, 1) / sqrt(sigma(0, 0) * sigma(1, 1));

  /* mvtnorm_C_mvtdst is defined in mvtnorm/inst/include/mvtnormAPI.h */
  mvtnorm_C_mvtdst(&n, &nu, lls.begin(), uls.begin(), infin, &corr, mu.begin(),
                   &maxpts, &abseps, &releps, &error, &value, &inform, &rnd);
  return value;
}


