# sppmix
Mixture models for spatial point pattern


This package implemented classes and methods for modeling spatial point process data using mixtures of normal component. 

Based on spatstat package's ppp class for point pattern and owin class for window, we implemented normmix class for dealing with 2D normal mixture.

Data augmentation MCMC (DAMCMC) and birth-death MCMC (BDMCMC) are the two main methods we have for estimating normal mixture to point pattern data.

The MCMC algorithms are implemented in C++ using Rcpp and RcppArmadillo, and it's significantly faster than some other implementations.

rgl package is used to create 3D intensity plots for intensity surface from normal mixture.

To learn more about sppmix, start with those vignettes: browseVignettes(package = "sppmix")

