#ifndef __SPPMIX_H__
#define __SPPMIX_H__

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdio.h>
#include <time.h>

using namespace arma;
using namespace Rcpp;

//Data augmentation, with truncation capabilities
//file: DAMCMC_sppmix.cpp
//FUNCION IS EXPOSED
//calls the functions: rMultinomial_sppmix,
//invmat2d_sppmix,rWishart_sppmix,ApproxMHRatiomu_sppmix,
//ApproxMHRatiosig_sppmix,ApproxCompMass_sppmix,
//rDirichlet_sppmix,densNormMixatx_sppmix
List DAMCMC2d_sppmix(mat const& data, 
                     vec const& xlims,
                     vec const& ylims,
                     int const& m,int const& L,
                     int const& burnin,int const& LL,
                     bool const& truncate);

//Simulation functions
//file: SimulFuncs_sppmix.cpp
ivec rMultinomial_sppmix(int const& n,vec const& ps);
List genNormMix_sppmix(int const& n,List const& mix);
double rUnifab_sppmix(double const& a,
                      double const& b);
double rUnif_sppmix();
mat rnorm2_sppmix(int n, vec mu,mat sigma);
mat rWishart_sppmix(int const& df, mat const& A);
int rDiscrete_sppmix(int const& start,
                     int const& end,
                     vec const& probs);
int rBinom_sppmix(int const& n,
                  double const& p);
double rGamma_sppmix(double const& a,
                     double const& b);
double rExp_sppmix(double const& a);
vec rDirichlet_sppmix(vec const& d);

//Helper functions
//file: HelperFuncs_sppmix.cpp
mat invmat2d_sppmix(mat const& A);
double densNormMixatx_sppmix(vec const& atx,
                             List const& mix);
mat dNormMix_sppmix(List const& mix, vec const& x,
                    vec const& y);
List GetStats_sppmix(vec const& gens,
                     double const& alpha);
mat ApproxAvgPostIntensity(List const& genmix,
                           int const& LL,
                           int const& burnin,
                           vec const& ticsx,
                           vec const& ticsy);
double ApproxCompMass_sppmix(int const& LL,vec const& ticsx,
                             vec const& ticsy,mat const& areas,
                             vec const& mu,
                             mat const& sig,mat const& siginv);
double ApproxMHRatiomu_sppmix(int const& LL,vec const& ticsx,
                              vec const& ticsy,mat const& areas,
                              vec const& curmu,vec const& propmu,
                              mat const& sig,mat const& siginv);
double ApproxMHRatiosig_sppmix(int const& LL,vec const& ticsx,
                               vec const& ticsy,mat const& areas,
                               vec const& mu1,mat const& propsigma,
                               mat const& sig,mat const& siginv);


#endif
