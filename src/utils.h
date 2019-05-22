#ifndef __BPGMM_UTILS__
#define __BPGMM_UTILS__

#include <RcppArmadillo.h>
// #include <iterator>
// #include <iostream>


arma::mat get_Z_mat(arma::vec ZOneDim, int m, int n);


arma::vec updatePostThetaY(int m, int n, Rcpp::S4 thetaYList);

Rcpp::NumericVector updatePostZ(int m,
                     int n,
                     Rcpp::S4 thetaYList);

Rcpp::List CalculateCxy(int m, int n, Rcpp::S4 hparam, Rcpp::S4 thetaYList,
                  arma::vec ZOneDim,
                  arma::vec qVec,
                  arma::mat X);

Rcpp::List CalculatePostLambdaPsy(Rcpp::S4 hparam,
                                  Rcpp::List CxyList,
                                  Rcpp::S4 thetaYList,
                                  arma::vec constraint);
#endif
