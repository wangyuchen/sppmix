// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/sppmix.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// ApproxAvgPostIntensity
mat ApproxAvgPostIntensity(List const& genmix, vec const& lamdas, int const& LL, int const& burnin, vec const& xlims, vec const& ylims);
RcppExport SEXP sppmix_ApproxAvgPostIntensity(SEXP genmixSEXP, SEXP lamdasSEXP, SEXP LLSEXP, SEXP burninSEXP, SEXP xlimsSEXP, SEXP ylimsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List const& >::type genmix(genmixSEXP);
    Rcpp::traits::input_parameter< vec const& >::type lamdas(lamdasSEXP);
    Rcpp::traits::input_parameter< int const& >::type LL(LLSEXP);
    Rcpp::traits::input_parameter< int const& >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< vec const& >::type xlims(xlimsSEXP);
    Rcpp::traits::input_parameter< vec const& >::type ylims(ylimsSEXP);
    __result = Rcpp::wrap(ApproxAvgPostIntensity(genmix, lamdas, LL, burnin, xlims, ylims));
    return __result;
END_RCPP
}
// ApproxCompMass_sppmix
double ApproxCompMass_sppmix(vec const& xlims, vec const& ylims, vec const& mu, mat const& sigma);
RcppExport SEXP sppmix_ApproxCompMass_sppmix(SEXP xlimsSEXP, SEXP ylimsSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< vec const& >::type xlims(xlimsSEXP);
    Rcpp::traits::input_parameter< vec const& >::type ylims(ylimsSEXP);
    Rcpp::traits::input_parameter< vec const& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< mat const& >::type sigma(sigmaSEXP);
    __result = Rcpp::wrap(ApproxCompMass_sppmix(xlims, ylims, mu, sigma));
    return __result;
END_RCPP
}
// ApproxMHRatiosig_sppmix
double ApproxMHRatiosig_sppmix(vec const& xlims, vec const& ylims, vec const& mu, mat const& cursigma, mat const& propsigma, int const& num);
RcppExport SEXP sppmix_ApproxMHRatiosig_sppmix(SEXP xlimsSEXP, SEXP ylimsSEXP, SEXP muSEXP, SEXP cursigmaSEXP, SEXP propsigmaSEXP, SEXP numSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< vec const& >::type xlims(xlimsSEXP);
    Rcpp::traits::input_parameter< vec const& >::type ylims(ylimsSEXP);
    Rcpp::traits::input_parameter< vec const& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< mat const& >::type cursigma(cursigmaSEXP);
    Rcpp::traits::input_parameter< mat const& >::type propsigma(propsigmaSEXP);
    Rcpp::traits::input_parameter< int const& >::type num(numSEXP);
    __result = Rcpp::wrap(ApproxMHRatiosig_sppmix(xlims, ylims, mu, cursigma, propsigma, num));
    return __result;
END_RCPP
}
// ApproxBayesianModelAvgIntensity_sppmix
mat ApproxBayesianModelAvgIntensity_sppmix(List const& genBDmix, vec const& lamdas, vec const& numcomp, vec const& distr_numcomp, int const& mincomp, int const& maxcomp, int const& LL, vec const& xlims, vec const& ylims);
RcppExport SEXP sppmix_ApproxBayesianModelAvgIntensity_sppmix(SEXP genBDmixSEXP, SEXP lamdasSEXP, SEXP numcompSEXP, SEXP distr_numcompSEXP, SEXP mincompSEXP, SEXP maxcompSEXP, SEXP LLSEXP, SEXP xlimsSEXP, SEXP ylimsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List const& >::type genBDmix(genBDmixSEXP);
    Rcpp::traits::input_parameter< vec const& >::type lamdas(lamdasSEXP);
    Rcpp::traits::input_parameter< vec const& >::type numcomp(numcompSEXP);
    Rcpp::traits::input_parameter< vec const& >::type distr_numcomp(distr_numcompSEXP);
    Rcpp::traits::input_parameter< int const& >::type mincomp(mincompSEXP);
    Rcpp::traits::input_parameter< int const& >::type maxcomp(maxcompSEXP);
    Rcpp::traits::input_parameter< int const& >::type LL(LLSEXP);
    Rcpp::traits::input_parameter< vec const& >::type xlims(xlimsSEXP);
    Rcpp::traits::input_parameter< vec const& >::type ylims(ylimsSEXP);
    __result = Rcpp::wrap(ApproxBayesianModelAvgIntensity_sppmix(genBDmix, lamdas, numcomp, distr_numcomp, mincomp, maxcomp, LL, xlims, ylims));
    return __result;
END_RCPP
}
// BDMCMC2d_sppmix
List BDMCMC2d_sppmix(int const& maxnumcomp, mat const& data, vec const& xlims, vec const& ylims, int const& L, int const& LL, bool const& truncate, double const& lamda, double const& lamdab, vec const& hyper);
RcppExport SEXP sppmix_BDMCMC2d_sppmix(SEXP maxnumcompSEXP, SEXP dataSEXP, SEXP xlimsSEXP, SEXP ylimsSEXP, SEXP LSEXP, SEXP LLSEXP, SEXP truncateSEXP, SEXP lamdaSEXP, SEXP lamdabSEXP, SEXP hyperSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int const& >::type maxnumcomp(maxnumcompSEXP);
    Rcpp::traits::input_parameter< mat const& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< vec const& >::type xlims(xlimsSEXP);
    Rcpp::traits::input_parameter< vec const& >::type ylims(ylimsSEXP);
    Rcpp::traits::input_parameter< int const& >::type L(LSEXP);
    Rcpp::traits::input_parameter< int const& >::type LL(LLSEXP);
    Rcpp::traits::input_parameter< bool const& >::type truncate(truncateSEXP);
    Rcpp::traits::input_parameter< double const& >::type lamda(lamdaSEXP);
    Rcpp::traits::input_parameter< double const& >::type lamdab(lamdabSEXP);
    Rcpp::traits::input_parameter< vec const& >::type hyper(hyperSEXP);
    __result = Rcpp::wrap(BDMCMC2d_sppmix(maxnumcomp, data, xlims, ylims, L, LL, truncate, lamda, lamdab, hyper));
    return __result;
END_RCPP
}
// DAMCMC2d_sppmix
List DAMCMC2d_sppmix(mat const& points, vec const& xlims, vec const& ylims, int const& m, int const& L, bool const& truncate, vec const& hyperparams);
RcppExport SEXP sppmix_DAMCMC2d_sppmix(SEXP pointsSEXP, SEXP xlimsSEXP, SEXP ylimsSEXP, SEXP mSEXP, SEXP LSEXP, SEXP truncateSEXP, SEXP hyperparamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< mat const& >::type points(pointsSEXP);
    Rcpp::traits::input_parameter< vec const& >::type xlims(xlimsSEXP);
    Rcpp::traits::input_parameter< vec const& >::type ylims(ylimsSEXP);
    Rcpp::traits::input_parameter< int const& >::type m(mSEXP);
    Rcpp::traits::input_parameter< int const& >::type L(LSEXP);
    Rcpp::traits::input_parameter< bool const& >::type truncate(truncateSEXP);
    Rcpp::traits::input_parameter< vec const& >::type hyperparams(hyperparamsSEXP);
    __result = Rcpp::wrap(DAMCMC2d_sppmix(points, xlims, ylims, m, L, truncate, hyperparams));
    return __result;
END_RCPP
}
// DAMCMC2dExtras_sppmix
List DAMCMC2dExtras_sppmix(mat const& points, vec const& xlims, vec const& ylims, int const& m, int const& L, int const& burnin, bool const& truncate, vec const& hyperparams);
RcppExport SEXP sppmix_DAMCMC2dExtras_sppmix(SEXP pointsSEXP, SEXP xlimsSEXP, SEXP ylimsSEXP, SEXP mSEXP, SEXP LSEXP, SEXP burninSEXP, SEXP truncateSEXP, SEXP hyperparamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< mat const& >::type points(pointsSEXP);
    Rcpp::traits::input_parameter< vec const& >::type xlims(xlimsSEXP);
    Rcpp::traits::input_parameter< vec const& >::type ylims(ylimsSEXP);
    Rcpp::traits::input_parameter< int const& >::type m(mSEXP);
    Rcpp::traits::input_parameter< int const& >::type L(LSEXP);
    Rcpp::traits::input_parameter< int const& >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< bool const& >::type truncate(truncateSEXP);
    Rcpp::traits::input_parameter< vec const& >::type hyperparams(hyperparamsSEXP);
    __result = Rcpp::wrap(DAMCMC2dExtras_sppmix(points, xlims, ylims, m, L, burnin, truncate, hyperparams));
    return __result;
END_RCPP
}
// Factorial_sppmix
double Factorial_sppmix(int x);
RcppExport SEXP sppmix_Factorial_sppmix(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    __result = Rcpp::wrap(Factorial_sppmix(x));
    return __result;
END_RCPP
}
// invmat2d_sppmix
mat invmat2d_sppmix(mat const& A);
RcppExport SEXP sppmix_invmat2d_sppmix(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< mat const& >::type A(ASEXP);
    __result = Rcpp::wrap(invmat2d_sppmix(A));
    return __result;
END_RCPP
}
// densNormMixatx_sppmix
double densNormMixatx_sppmix(vec const& atx, List const& mix);
RcppExport SEXP sppmix_densNormMixatx_sppmix(SEXP atxSEXP, SEXP mixSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< vec const& >::type atx(atxSEXP);
    Rcpp::traits::input_parameter< List const& >::type mix(mixSEXP);
    __result = Rcpp::wrap(densNormMixatx_sppmix(atx, mix));
    return __result;
END_RCPP
}
// densNormMix_atxy_sppmix
vec densNormMix_atxy_sppmix(mat const& atxy, List const& mix);
RcppExport SEXP sppmix_densNormMix_atxy_sppmix(SEXP atxySEXP, SEXP mixSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< mat const& >::type atxy(atxySEXP);
    Rcpp::traits::input_parameter< List const& >::type mix(mixSEXP);
    __result = Rcpp::wrap(densNormMix_atxy_sppmix(atxy, mix));
    return __result;
END_RCPP
}
// dNormMix_sppmix
mat dNormMix_sppmix(List const& mix, vec const& x, vec const& y);
RcppExport SEXP sppmix_dNormMix_sppmix(SEXP mixSEXP, SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List const& >::type mix(mixSEXP);
    Rcpp::traits::input_parameter< vec const& >::type x(xSEXP);
    Rcpp::traits::input_parameter< vec const& >::type y(ySEXP);
    __result = Rcpp::wrap(dNormMix_sppmix(mix, x, y));
    return __result;
END_RCPP
}
// Permute_vec_sppmix
vec Permute_vec_sppmix(vec const& oldvec, vec const& perm);
RcppExport SEXP sppmix_Permute_vec_sppmix(SEXP oldvecSEXP, SEXP permSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< vec const& >::type oldvec(oldvecSEXP);
    Rcpp::traits::input_parameter< vec const& >::type perm(permSEXP);
    __result = Rcpp::wrap(Permute_vec_sppmix(oldvec, perm));
    return __result;
END_RCPP
}
// Permute_mat_sppmix
mat Permute_mat_sppmix(mat const& oldmat, vec const& perm);
RcppExport SEXP sppmix_Permute_mat_sppmix(SEXP oldmatSEXP, SEXP permSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< mat const& >::type oldmat(oldmatSEXP);
    Rcpp::traits::input_parameter< vec const& >::type perm(permSEXP);
    __result = Rcpp::wrap(Permute_mat_sppmix(oldmat, perm));
    return __result;
END_RCPP
}
// GetAllPermutations_sppmix
mat GetAllPermutations_sppmix(int const& m);
RcppExport SEXP sppmix_GetAllPermutations_sppmix(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int const& >::type m(mSEXP);
    __result = Rcpp::wrap(GetAllPermutations_sppmix(m));
    return __result;
END_RCPP
}
// GetAPermutation_sppmix
vec GetAPermutation_sppmix(int const& m, int const& which);
RcppExport SEXP sppmix_GetAPermutation_sppmix(SEXP mSEXP, SEXP whichSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int const& >::type m(mSEXP);
    Rcpp::traits::input_parameter< int const& >::type which(whichSEXP);
    __result = Rcpp::wrap(GetAPermutation_sppmix(m, which));
    return __result;
END_RCPP
}
// GetGrid_sppmix
List GetGrid_sppmix(int const& len, vec const& xlims, vec const& ylims);
RcppExport SEXP sppmix_GetGrid_sppmix(SEXP lenSEXP, SEXP xlimsSEXP, SEXP ylimsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int const& >::type len(lenSEXP);
    Rcpp::traits::input_parameter< vec const& >::type xlims(xlimsSEXP);
    Rcpp::traits::input_parameter< vec const& >::type ylims(ylimsSEXP);
    __result = Rcpp::wrap(GetGrid_sppmix(len, xlims, ylims));
    return __result;
END_RCPP
}
// EqVec_sppmix
bool EqVec_sppmix(vec const& v1, vec const& v2, double const& tol);
RcppExport SEXP sppmix_EqVec_sppmix(SEXP v1SEXP, SEXP v2SEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< vec const& >::type v1(v1SEXP);
    Rcpp::traits::input_parameter< vec const& >::type v2(v2SEXP);
    Rcpp::traits::input_parameter< double const& >::type tol(tolSEXP);
    __result = Rcpp::wrap(EqVec_sppmix(v1, v2, tol));
    return __result;
END_RCPP
}
// logGammaFunc_sppmix
double logGammaFunc_sppmix(double const& x);
RcppExport SEXP sppmix_logGammaFunc_sppmix(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double const& >::type x(xSEXP);
    __result = Rcpp::wrap(logGammaFunc_sppmix(x));
    return __result;
END_RCPP
}
// GammaFunc_sppmix
double GammaFunc_sppmix(double const& x);
RcppExport SEXP sppmix_GammaFunc_sppmix(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double const& >::type x(xSEXP);
    __result = Rcpp::wrap(GammaFunc_sppmix(x));
    return __result;
END_RCPP
}
// dDirichlet_sppmix
double dDirichlet_sppmix(vec const& ps, vec const& ds);
RcppExport SEXP sppmix_dDirichlet_sppmix(SEXP psSEXP, SEXP dsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< vec const& >::type ps(psSEXP);
    Rcpp::traits::input_parameter< vec const& >::type ds(dsSEXP);
    __result = Rcpp::wrap(dDirichlet_sppmix(ps, ds));
    return __result;
END_RCPP
}
// SumVec_sppmix
double SumVec_sppmix(vec const& v, int const& start, int const& end);
RcppExport SEXP sppmix_SumVec_sppmix(SEXP vSEXP, SEXP startSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< vec const& >::type v(vSEXP);
    Rcpp::traits::input_parameter< int const& >::type start(startSEXP);
    Rcpp::traits::input_parameter< int const& >::type end(endSEXP);
    __result = Rcpp::wrap(SumVec_sppmix(v, start, end));
    return __result;
END_RCPP
}
// SubstituteVec_sppmix
vec SubstituteVec_sppmix(vec v, vec const& subv, int const& start);
RcppExport SEXP sppmix_SubstituteVec_sppmix(SEXP vSEXP, SEXP subvSEXP, SEXP startSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< vec >::type v(vSEXP);
    Rcpp::traits::input_parameter< vec const& >::type subv(subvSEXP);
    Rcpp::traits::input_parameter< int const& >::type start(startSEXP);
    __result = Rcpp::wrap(SubstituteVec_sppmix(v, subv, start));
    return __result;
END_RCPP
}
// SubVec_sppmix
vec SubVec_sppmix(vec const& v, int const& start, int const& end);
RcppExport SEXP sppmix_SubVec_sppmix(SEXP vSEXP, SEXP startSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< vec const& >::type v(vSEXP);
    Rcpp::traits::input_parameter< int const& >::type start(startSEXP);
    Rcpp::traits::input_parameter< int const& >::type end(endSEXP);
    __result = Rcpp::wrap(SubVec_sppmix(v, start, end));
    return __result;
END_RCPP
}
// GetMixtureMaxz_sppmix
double GetMixtureMaxz_sppmix(List const& genmix, int const& len, vec const& xlims, vec const& ylims);
RcppExport SEXP sppmix_GetMixtureMaxz_sppmix(SEXP genmixSEXP, SEXP lenSEXP, SEXP xlimsSEXP, SEXP ylimsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List const& >::type genmix(genmixSEXP);
    Rcpp::traits::input_parameter< int const& >::type len(lenSEXP);
    Rcpp::traits::input_parameter< vec const& >::type xlims(xlimsSEXP);
    Rcpp::traits::input_parameter< vec const& >::type ylims(ylimsSEXP);
    __result = Rcpp::wrap(GetMixtureMaxz_sppmix(genmix, len, xlims, ylims));
    return __result;
END_RCPP
}
// MakeMixtureList_sppmix
List MakeMixtureList_sppmix(List const& gens_list, int const& burnin);
RcppExport SEXP sppmix_MakeMixtureList_sppmix(SEXP gens_listSEXP, SEXP burninSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List const& >::type gens_list(gens_listSEXP);
    Rcpp::traits::input_parameter< int const& >::type burnin(burninSEXP);
    __result = Rcpp::wrap(MakeMixtureList_sppmix(gens_list, burnin));
    return __result;
END_RCPP
}
// CheckInWindow_sppmix
List CheckInWindow_sppmix(mat const& points, vec const& xlims, vec const& ylims, bool const& truncate);
RcppExport SEXP sppmix_CheckInWindow_sppmix(SEXP pointsSEXP, SEXP xlimsSEXP, SEXP ylimsSEXP, SEXP truncateSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< mat const& >::type points(pointsSEXP);
    Rcpp::traits::input_parameter< vec const& >::type xlims(xlimsSEXP);
    Rcpp::traits::input_parameter< vec const& >::type ylims(ylimsSEXP);
    Rcpp::traits::input_parameter< bool const& >::type truncate(truncateSEXP);
    __result = Rcpp::wrap(CheckInWindow_sppmix(points, xlims, ylims, truncate));
    return __result;
END_RCPP
}
// GetMax_sppmix
List GetMax_sppmix(vec const& v);
RcppExport SEXP sppmix_GetMax_sppmix(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< vec const& >::type v(vSEXP);
    __result = Rcpp::wrap(GetMax_sppmix(v));
    return __result;
END_RCPP
}
// GetStats_sppmix
List GetStats_sppmix(vec const& gens, double const& alpha);
RcppExport SEXP sppmix_GetStats_sppmix(SEXP gensSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< vec const& >::type gens(gensSEXP);
    Rcpp::traits::input_parameter< double const& >::type alpha(alphaSEXP);
    __result = Rcpp::wrap(GetStats_sppmix(gens, alpha));
    return __result;
END_RCPP
}
// GetRealiz_ps_sppmix
vec GetRealiz_ps_sppmix(List const& allgens, int const& realiz);
RcppExport SEXP sppmix_GetRealiz_ps_sppmix(SEXP allgensSEXP, SEXP realizSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List const& >::type allgens(allgensSEXP);
    Rcpp::traits::input_parameter< int const& >::type realiz(realizSEXP);
    __result = Rcpp::wrap(GetRealiz_ps_sppmix(allgens, realiz));
    return __result;
END_RCPP
}
// GetRealiz_mus_sppmix
mat GetRealiz_mus_sppmix(List const& allgens, int const& realiz);
RcppExport SEXP sppmix_GetRealiz_mus_sppmix(SEXP allgensSEXP, SEXP realizSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List const& >::type allgens(allgensSEXP);
    Rcpp::traits::input_parameter< int const& >::type realiz(realizSEXP);
    __result = Rcpp::wrap(GetRealiz_mus_sppmix(allgens, realiz));
    return __result;
END_RCPP
}
// GetRealiz_sigmas_sppmix
mat GetRealiz_sigmas_sppmix(List const& allgens, int const& realiz);
RcppExport SEXP sppmix_GetRealiz_sigmas_sppmix(SEXP allgensSEXP, SEXP realizSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List const& >::type allgens(allgensSEXP);
    Rcpp::traits::input_parameter< int const& >::type realiz(realizSEXP);
    __result = Rcpp::wrap(GetRealiz_sigmas_sppmix(allgens, realiz));
    return __result;
END_RCPP
}
// PostGenGetBestPerm_sppmix
List PostGenGetBestPerm_sppmix(List const& allgens);
RcppExport SEXP sppmix_PostGenGetBestPerm_sppmix(SEXP allgensSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List const& >::type allgens(allgensSEXP);
    __result = Rcpp::wrap(PostGenGetBestPerm_sppmix(allgens));
    return __result;
END_RCPP
}
// GetAllMeans_sppmix
List GetAllMeans_sppmix(List const& allgens, int const& burnin);
RcppExport SEXP sppmix_GetAllMeans_sppmix(SEXP allgensSEXP, SEXP burninSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List const& >::type allgens(allgensSEXP);
    Rcpp::traits::input_parameter< int const& >::type burnin(burninSEXP);
    __result = Rcpp::wrap(GetAllMeans_sppmix(allgens, burnin));
    return __result;
END_RCPP
}
// GetCompDistr_sppmix
vec GetCompDistr_sppmix(vec const& numcomp, int const& maxnumcomp);
RcppExport SEXP sppmix_GetCompDistr_sppmix(SEXP numcompSEXP, SEXP maxnumcompSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< vec const& >::type numcomp(numcompSEXP);
    Rcpp::traits::input_parameter< int const& >::type maxnumcomp(maxnumcompSEXP);
    __result = Rcpp::wrap(GetCompDistr_sppmix(numcomp, maxnumcomp));
    return __result;
END_RCPP
}
// GetBDCompRealiz_sppmix
List GetBDCompRealiz_sppmix(List const& genBDmix, vec const& genlamdas, vec const& numcomp, int const& comp);
RcppExport SEXP sppmix_GetBDCompRealiz_sppmix(SEXP genBDmixSEXP, SEXP genlamdasSEXP, SEXP numcompSEXP, SEXP compSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List const& >::type genBDmix(genBDmixSEXP);
    Rcpp::traits::input_parameter< vec const& >::type genlamdas(genlamdasSEXP);
    Rcpp::traits::input_parameter< vec const& >::type numcomp(numcompSEXP);
    Rcpp::traits::input_parameter< int const& >::type comp(compSEXP);
    __result = Rcpp::wrap(GetBDCompRealiz_sppmix(genBDmix, genlamdas, numcomp, comp));
    return __result;
END_RCPP
}
// GetAvgLabelsDiscrete2Multinomial_sppmix
mat GetAvgLabelsDiscrete2Multinomial_sppmix(mat       const& genzs, int const& m);
RcppExport SEXP sppmix_GetAvgLabelsDiscrete2Multinomial_sppmix(SEXP genzsSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< mat       const& >::type genzs(genzsSEXP);
    Rcpp::traits::input_parameter< int const& >::type m(mSEXP);
    __result = Rcpp::wrap(GetAvgLabelsDiscrete2Multinomial_sppmix(genzs, m));
    return __result;
END_RCPP
}
// Check4LabelSwitching_sppmix
bool Check4LabelSwitching_sppmix(vec const& chain);
RcppExport SEXP sppmix_Check4LabelSwitching_sppmix(SEXP chainSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< vec const& >::type chain(chainSEXP);
    __result = Rcpp::wrap(Check4LabelSwitching_sppmix(chain));
    return __result;
END_RCPP
}
// ApproxBivNormProb_sppmix
double ApproxBivNormProb_sppmix(vec const& xlims, vec const& ylims, vec const& mu, mat const& sigma, int type);
RcppExport SEXP sppmix_ApproxBivNormProb_sppmix(SEXP xlimsSEXP, SEXP ylimsSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< vec const& >::type xlims(xlimsSEXP);
    Rcpp::traits::input_parameter< vec const& >::type ylims(ylimsSEXP);
    Rcpp::traits::input_parameter< vec const& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< mat const& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    __result = Rcpp::wrap(ApproxBivNormProb_sppmix(xlims, ylims, mu, sigma, type));
    return __result;
END_RCPP
}
// rnorm2_sppmix
mat rnorm2_sppmix(int const& n, vec const& mu, mat const& sigma);
RcppExport SEXP sppmix_rnorm2_sppmix(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int const& >::type n(nSEXP);
    Rcpp::traits::input_parameter< vec const& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< mat const& >::type sigma(sigmaSEXP);
    __result = Rcpp::wrap(rnorm2_sppmix(n, mu, sigma));
    return __result;
END_RCPP
}
// rWishart_sppmix
mat rWishart_sppmix(int const& df, mat const& A);
RcppExport SEXP sppmix_rWishart_sppmix(SEXP dfSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int const& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< mat const& >::type A(ASEXP);
    __result = Rcpp::wrap(rWishart_sppmix(df, A));
    return __result;
END_RCPP
}
// rDiscrete_sppmix
int rDiscrete_sppmix(int const& start, vec const& probs);
RcppExport SEXP sppmix_rDiscrete_sppmix(SEXP startSEXP, SEXP probsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int const& >::type start(startSEXP);
    Rcpp::traits::input_parameter< vec const& >::type probs(probsSEXP);
    __result = Rcpp::wrap(rDiscrete_sppmix(start, probs));
    return __result;
END_RCPP
}
// rBinom_sppmix
int rBinom_sppmix(int const& n, double const& p);
RcppExport SEXP sppmix_rBinom_sppmix(SEXP nSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int const& >::type n(nSEXP);
    Rcpp::traits::input_parameter< double const& >::type p(pSEXP);
    __result = Rcpp::wrap(rBinom_sppmix(n, p));
    return __result;
END_RCPP
}
// rExp_sppmix
double rExp_sppmix(double const& a);
RcppExport SEXP sppmix_rExp_sppmix(SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double const& >::type a(aSEXP);
    __result = Rcpp::wrap(rExp_sppmix(a));
    return __result;
END_RCPP
}
// rDirichlet_sppmix
vec rDirichlet_sppmix(vec const& d);
RcppExport SEXP sppmix_rDirichlet_sppmix(SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< vec const& >::type d(dSEXP);
    __result = Rcpp::wrap(rDirichlet_sppmix(d));
    return __result;
END_RCPP
}
// rMultinomial_sppmix
ivec rMultinomial_sppmix(int const& n, vec const& ps);
RcppExport SEXP sppmix_rMultinomial_sppmix(SEXP nSEXP, SEXP psSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int const& >::type n(nSEXP);
    Rcpp::traits::input_parameter< vec const& >::type ps(psSEXP);
    __result = Rcpp::wrap(rMultinomial_sppmix(n, ps));
    return __result;
END_RCPP
}
// rNormMix_sppmix
List rNormMix_sppmix(int const& lamda, List const& mix);
RcppExport SEXP sppmix_rNormMix_sppmix(SEXP lamdaSEXP, SEXP mixSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int const& >::type lamda(lamdaSEXP);
    Rcpp::traits::input_parameter< List const& >::type mix(mixSEXP);
    __result = Rcpp::wrap(rNormMix_sppmix(lamda, mix));
    return __result;
END_RCPP
}
// rPerm_sppmix
vec rPerm_sppmix(int const& n);
RcppExport SEXP sppmix_rPerm_sppmix(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int const& >::type n(nSEXP);
    __result = Rcpp::wrap(rPerm_sppmix(n));
    return __result;
END_RCPP
}
