// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <iostream>
#include "utils.h"
using namespace Rcpp;


// [[Rcpp::export]]
arma::vec update_PostThetaY(int m, int n, Rcpp::S4 thetaYList) {

  arma::vec ZOneDim(n);
  return ZOneDim;
}
