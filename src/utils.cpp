// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <cmath>
#include "utils.h"
using namespace Rcpp;


// [[Rcpp::export]]
arma::mat get_Z_mat(arma::vec ZOneDim, int m, int n){

  if (m < 1) {
    Rcpp::stop("m must be a positive integer");
  }
  if (n < 1) {
    Rcpp::stop("n must be a positive integer");
  }
  if (ZOneDim.n_elem != static_cast<arma::uword>(n)) {
    Rcpp::stop("length of ZOneDim must equal n");
  }

  arma::mat Zmat = arma::zeros<arma::mat>(m, n);

  for (int j = 0; j < n; j++) {
    double label = ZOneDim(j);
    if (!std::isfinite(label) || label < 1 || label > m || label != std::floor(label)) {
      Rcpp::stop("cluster labels must be integers in 1:m");
    }
    Zmat(static_cast<arma::uword>(label - 1), j) = 1;
  }
  return(Zmat);
}

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
arma::vec dmvnrm_arma(arma::mat x,
                      arma::rowvec mean,
                      arma::mat sigma,
                      bool logd) {
  int n = x.n_rows;
  int xdim = x.n_cols;

  if (xdim < 1) {
    Rcpp::stop("x must have at least one column");
  }
  if (mean.n_elem != static_cast<arma::uword>(xdim)) {
    Rcpp::stop("mean length must equal the number of columns in x");
  }
  if (sigma.n_rows != sigma.n_cols || sigma.n_rows != static_cast<arma::uword>(xdim)) {
    Rcpp::stop("sigma must be a square matrix with dimension matching x");
  }
  if (!x.is_finite() || !mean.is_finite() || !sigma.is_finite()) {
    Rcpp::stop("x, mean, and sigma must contain only finite values");
  }

  arma::vec out(n);
  arma::mat sigma_chol;
  if (!arma::chol(sigma_chol, sigma)) {
    Rcpp::stop("sigma must be positive definite");
  }

  arma::mat rooti = arma::trans(arma::inv(trimatu(sigma_chol)));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;

  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;
  }

  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

// [[Rcpp::export]]
double calculate_Ratio(double logDeno, arma::vec logNume){

  if (!std::isfinite(logDeno)) {
    Rcpp::stop("logDeno must be finite");
  }
  if (logNume.n_elem < 1) {
    Rcpp::stop("logNume must contain at least one value");
  }
  if (!logNume.is_finite()) {
    Rcpp::stop("logNume must contain only finite values");
  }

  int n = logNume.n_elem;
  double maxNume = arma::max(logNume);
  double transDeno = logDeno - maxNume;

  arma::vec repMaxNume = rep(maxNume,n);

  arma::vec transNume = logNume - repMaxNume;
  double ratio = exp(transDeno)/sum(exp(transNume));

  return(ratio);
}


