#ifndef __BPGMM_UTILS__
#define __BPGMM_UTILS__

#include <RcppArmadillo.h>
// #include <iterator>
// #include <iostream>

arma::vec dmvnrm_arma(arma::mat x,
                      arma::rowvec mean,
                      arma::mat sigma,
                      bool logd);


arma::mat get_Z_mat(arma::vec ZOneDim, int m, int n);

double calculate_Ratio(double logDeno, arma::vec logNume);

arma::vec update_PostThetaY(int m, int n, Rcpp::S4 thetaYList);

Rcpp::NumericVector updatePost_Z( arma::mat X,int m,
                     int n,
                     Rcpp::S4 thetaYList);

Rcpp::List CalculateCxy(int m, int n, Rcpp::S4 hparam, Rcpp::S4 thetaYList,
                  arma::vec ZOneDim,
                  arma::vec qVec,
                  arma::mat X);

Rcpp::List Calculate_PostLambdaPsy(Rcpp::S4 hparam,
                                   Rcpp::List CxyList,
                                   Rcpp::S4 thetaYList,
                                   arma::vec constraint);

// S4 update_Hyperparameter(
//     int m,
//     int p,
//     int q,
//     Rcpp::S4 hparam,
//     Rcpp::S4 thetaYList,
//     arma::vec dVec,
//     arma::vec sVec
// );
#endif
