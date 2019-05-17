// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include "utils.h"
using namespace Rcpp;


// [[Rcpp::export]]
arma::mat get_Z_mat(arma::vec ZOneDim, int m, int n){

  arma::mat Zmat = arma::zeros<arma::mat>(m, n);

  for (int j = 0; j < n; j++) {
    Zmat(ZOneDim(j) - 1, j) = 1;
  }
  return(Zmat);
}
