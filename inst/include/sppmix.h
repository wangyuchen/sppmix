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
double rUnif_sppmix();
double rUnifab_sppmix(double const& a,
                      double const& b);
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
ivec rMultinomial_sppmix(int const& n,vec const& ps);
List rNormMix_sppmix(int const& lamda,List const& mix);
vec rPerm_sppmix(int const& n);

//Approximation functions
//file: ApproxFuncs_sppmix.cpp
mat ApproxAvgPostIntensity(List const& genmix,
    vec const& lamdas,int const& LL,int const& burnin,
    vec const& ticsx,vec const& ticsy);
double ApproxCompMass_sppmix(int const& LL,vec const& ticsx,
    vec const& ticsy,mat const& areas,vec const& mu,
    mat const& sig,mat const& siginv);
double ApproxMHRatiomu_sppmix(int const& LL,vec const& ticsx,
    vec const& ticsy,mat const& areas,
    vec const& curmu,vec const& propmu,
    mat const& sig,mat const& siginv);
double ApproxMHRatiosig_sppmix(int const& LL,vec const& ticsx,
    vec const& ticsy,mat const& areas,
    vec const& mu1,mat const& propsigma,
    mat const& sig,mat const& siginv);

//Operations on posterior realizations
//file: OpsPostGens_sppmix.cpp
List GetStats_sppmix(vec const& gens,double const& alpha);
vec GetRealiz_ps_sppmix(List const& allgens,
                        int const& realiz);
mat GetRealiz_mus_sppmix(List const& allgens,
                         int const& realiz);
mat GetRealiz_sigmas_sppmix(List const& allgens,
                            int const& realiz);
List PostGenGetBestPerm_sppmix(List const& allgens);
List GetAllMeans_sppmix(List const& allgens,int const& burnin);

//Helper functions
//file: HelperFuncs_sppmix.cpp
double Factorial_sppmix(int x);
mat invmat2d_sppmix(mat const& A);
double densNormMixatx_sppmix(vec const& atx,List const& mix);
mat dNormMix_sppmix(List const& mix, vec const& x,vec const& y);
vec Permute_vec_sppmix(vec const& oldvec,vec const& perm);
mat Permute_mat_sppmix(mat const& oldmat,vec const& perm);
mat GetAllPermutations_sppmix(int const& m);
vec GetAPermutation_sppmix(int const& m,int const& which);
List GetGrid_sppmix(int const& len,vec const& mins,vec const& maxs);
bool EqVec_sppmix(vec const& v1,vec const& v2,double const& tol=0.000001);
double dDirichlet_sppmix(vec const& ps,vec const& ds);
double logGammaFunc_sppmix(double const& x);
double GammaFunc_sppmix(double const& x);

#endif
