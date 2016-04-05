#ifndef __sppmix_h__
#define __sppmix_h__

// It's ok to use #include in headers, but the header guard is necessary.


// Don't use namespace in headers, only in .cpp files.
// Always use explicit naming in headers.

arma::vec dbvnorm(arma::mat x, arma::rowvec mean,
                  arma::mat sigma, bool logd = false);

#endif
