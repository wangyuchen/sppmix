#include <Rcpp.h>
#include <mvtnormAPI.h>
#include "sppmix.h"

using namespace Rcpp;

// [[Rcpp::export]]
double pbvnorm(NumericVector const& xlims, NumericVector const& ylims,
               NumericVector const& mu, NumericMatrix const& sigma) {
  NumericVector lls(2), uls(2);
  lls(0) = (xlims(0) - mu(0)) / sqrt(sigma(0, 0));
  lls(1) = (ylims(0) - mu(1)) / sqrt(sigma(1, 1));

  uls(0) = (xlims(1) - mu(0)) / sqrt(sigma(0, 0));
  uls(1) = (ylims(1) - mu(1)) / sqrt(sigma(1, 1));

  NumericVector muxy = NumericVector::create(0, 0);

  int n = 2, nu = 0, maxpts = 2000, inform;

  //INFIN INTEGER, array of integration limits flags:
  //if INFIN(I) < 0, Ith limits are (-infinity, infinity);
  //if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
  //if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
  //if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
  IntegerVector infin = IntegerVector::create(2, 2);

  double abseps = 1/1000, releps = 1 / 1000, error, value;
  int rnd = 0;
  double corr = sigma(0, 1) / sqrt(sigma(0, 0) * sigma(1, 1));

  /* mvtnorm_C_mvtdst is defined in mvtnorm/inst/include/mvtnormAPI.h */
  mvtnorm_C_mvtdst(&n, &nu, lls.begin(), uls.begin(),
                   infin.begin(), &corr, muxy.begin(),
                   &maxpts, &abseps, &releps, &error,
                   &value, &inform, &rnd);
                   return value;
}



//' Approximate density of a normal mixture over a 2d domain.
//'
//' Approximate the density of each component in a normal mixture within the
//' domain using multivariate normal density function.
//'
//' @param mix An object of class \code{normmix}
//' @param xlim,ylim Vector of length two. Mixture density are estimated within
//' this range.
//' @importFrom mvtnorm dmvnorm
//'
//' @return A numerical vector corresponding to the density of each component
//'  within the window.
//' @export
// [[Rcpp::export]]
NumericVector approx_normmix(List mix, NumericVector xlim,
                             NumericVector ylim) {
  if (!mix.inherits("normmix"))
    stop("mix must be of class normmix or intensity surface");

  int m = as<int>(mix["m"]);

  List mus = as<List>(mix["mus"]);
  List sigmas = as<List>(mix["sigmas"]);

  NumericVector approx(m);
  NumericVector mu(2);
  NumericMatrix sigma(2, 2);

  for (int i = 0; i < m; ++i) {
    mu = as<NumericVector>(mus[i]);
    sigma = as<NumericMatrix>(sigmas[i]);
    approx(i) = pbvnorm(xlim, ylim, mu, sigma);
  }

  return approx;
}











