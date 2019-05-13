// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
arma::mat get_Z_mat(arma::vec ZOneDim, int m, int n){

  arma::mat Zmat = arma::zeros<arma::mat>(m, n);

  for (int j = 0; j < n; j++) {
    Zmat(ZOneDim(j) - 1, j) = 1;
  }
  return(Zmat);
}



/*** R
ZOneDim1 <- c(1,2,3,1,2,2)
ZOneDim2 <- c(1,1,3,1,2,2)

Zmat1 <- matrix(c(1,0,0,0,1,0,0,0,1,1,0,0,0,1,0,0,1,0),3,6)
Zmat2 <- matrix(c(1,0,0,1,0,0,0,0,1,1,0,0,0,1,0,0,1,0),3,6)

testthat::expect_equal(get_Z_mat(ZOneDim1,3,6),Zmat1)
testthat::expect_equal(get_Z_mat(ZOneDim2,3,6),Zmat2)
*/
