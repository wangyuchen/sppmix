#include "sppmix.h"
#include <mvtnormAPI.h>


// [[Rcpp::export]]
double pbvnorm(Rcpp::NumericVector const& xlims,
               Rcpp::NumericVector const& ylims,
               Rcpp::NumericVector const& mu,
               Rcpp::NumericMatrix const& sigma) {
  Rcpp::NumericVector lls(2), uls(2);
  lls(0) = (xlims(0) - mu(0)) / sqrt(sigma(0, 0));
  lls(1) = (ylims(0) - mu(1)) / sqrt(sigma(1, 1));

  uls(0) = (xlims(1) - mu(0)) / sqrt(sigma(0, 0));
  uls(1) = (ylims(1) - mu(1)) / sqrt(sigma(1, 1));

  Rcpp::NumericVector muxy = Rcpp::NumericVector::create(0, 0);

  int n = 2, nu = 0, maxpts = 2000, inform;

  //INFIN INTEGER, array of integration limits flags:
  //if INFIN(I) < 0, Ith limits are (-infinity, infinity);
  //if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
  //if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
  //if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
  Rcpp::IntegerVector infin = Rcpp::IntegerVector::create(2, 2);

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
//' @return A numerical vector corresponding to the density of each
//'  component
//'  within the window.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector approx_normmix(Rcpp::List mix, Rcpp::NumericVector xlim,
                                   Rcpp::NumericVector ylim) {
  if (!mix.inherits("normmix"))
    Rcpp::stop("mix must be of class normmix or intensity surface");

  int m = Rcpp::as<int>(mix["m"]);

  Rcpp::List mus = Rcpp::as<Rcpp::List>(mix["mus"]);
  Rcpp::List sigmas = Rcpp::as<Rcpp::List>(mix["sigmas"]);

  Rcpp::NumericVector approx(m);
  Rcpp::NumericVector mu(2);
  Rcpp::NumericMatrix sigma(2, 2);

  for (int i = 0; i < m; ++i) {
    mu = Rcpp::as<Rcpp::NumericVector>(mus[i]);
    sigma = Rcpp::as<Rcpp::NumericMatrix>(sigmas[i]);
    approx(i) = pbvnorm(xlim, ylim, mu, sigma);
  }

  return approx;
}


//' Calculate density of normal mixture.
//'
//' When a \code{normmix} object is given, this function estimates density
//' values within the domain. When a \code{intensity_surface} is given, it
//' multiplies the estimated density values with intensity and returns the
//' estimated intensity.
//'
//' @param mix An object of class \code{normmix} or \code{intensity_surface}.
//' When an intensity surface is given, it estimates the intensity instead of
//' density within the domain.
//' @param xlim,ylim The range within which the density or intensity is
//' calculated.
//' @param L Length of grid on each axis. The density is calculated on a L * L
//' grid.
//' @param truncate Whether to truncate the density to be within the domain.
//'
//' @return An matrix with the first two columns the grid coordinates and the
//' third column the density (intensity) value.
//'
//'
//' @examples
//' est_density <- dnormmix(demo_mix)
//' est_density <- dnormmix(demo_intsurf)
//' @export
// [[Rcpp::export]]
arma::mat dnormmix(Rcpp::List mix,
                   Rcpp::NumericVector xlim = Rcpp::NumericVector::create(0, 1),
                   Rcpp::NumericVector ylim = Rcpp::NumericVector::create(0, 1),
                   int L = 128, bool truncate = true) {
  if (!mix.inherits("normmix"))
    Rcpp::stop("mix must be of class normmix or intensity surface");

  int m = Rcpp::as<int>(mix["m"]);
  arma::vec lps = log(Rcpp::as<arma::vec>(mix["ps"]));

  Rcpp::List mus = Rcpp::as<Rcpp::List>(mix["mus"]);
  Rcpp::List sigmas = Rcpp::as<Rcpp::List>(mix["sigmas"]);

  if (mix.inherits("intensity_surface")) {
    xlim = Rcpp::as<Rcpp::NumericVector>(
      Rcpp::as<Rcpp::List>(mix["window"])["xrange"]);
    ylim = Rcpp::as<Rcpp::NumericVector>(
      Rcpp::as<Rcpp::List>(mix["window"])["yrange"]);

    double intensity = Rcpp::as<double>(mix["intensity"]);
    lps = lps + log(intensity); // calculate intensity instead of density
  }

  arma::mat locs(L * L, 2);
  for (int i = 0; i < L; ++i) {
    for (int j = 0; j < L; ++j) {
      locs(i * L + j, 0) = xlim[0] + i * (xlim[1] - xlim[0]) / (L - 1);
      locs(i * L + j, 1) = ylim[0] + j * (ylim[1] - ylim[0]) / (L - 1);
    }
  }

  arma::mat lden(L * L, m);
  for (int k = 0; k < m; ++k) {
    lden.col(k) = dbvnorm(locs, Rcpp::as<arma::rowvec>(mus[k]),
             Rcpp::as<arma::mat>(sigmas[k]), true);
  }

  if (truncate) {
    arma::vec approx = Rcpp::as<arma::vec>(approx_normmix(mix, xlim, ylim));
    lps = lps - log(approx);
  }

  arma::mat ret = join_rows(locs, exp(lden) * exp(lps));
  return ret;
}


