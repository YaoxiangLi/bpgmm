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
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
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

  int n = logNume.n_elem;
  double maxNume = arma::max(logNume);
  double transDeno = logDeno - maxNume;

  arma::vec repMaxNume = rep(maxNume,n);

  arma::vec transNume = logNume - repMaxNume;
  double ratio = exp(transDeno)/sum(exp(transNume));

  return(ratio);
}



