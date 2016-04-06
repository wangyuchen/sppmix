#ifndef __sppmix_h__
#define __sppmix_h__

#include <RcppArmadillo.h>

arma::vec dbvnorm(arma::mat x, arma::rowvec mean,
                  arma::mat sigma, bool logd);

#endif
