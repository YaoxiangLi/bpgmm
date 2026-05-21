// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include "utils.h"
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::List Calculate_Cxy(int m,
                         int n,
                         Rcpp::S4 hparam,
                         Rcpp::S4 thetaYList,
                         arma::vec ZOneDim,
                         arma::vec qVec,
                         arma::mat X){

  validate_positive_int(m, "m");
  validate_positive_int(n, "n");
  if (X.n_rows < 1) {
    Rcpp::stop("X must have at least one row");
  }
  if (X.n_cols != static_cast<arma::uword>(n)) {
    Rcpp::stop("n must equal the number of columns in X");
  }
  validate_finite_matrix(X, "X");
  validate_q_vec(qVec, m);

  double alpha1 = get_positive_finite_slot(hparam, "alpha1");
  double alpha2 = get_positive_finite_slot(hparam, "alpha2");

  List Y      = thetaYList.slot("Y");
  List lambda = thetaYList.slot("lambda");
  List M      = thetaYList.slot("M");
  List psy    = thetaYList.slot("psy");

  if (Y.size() < m || lambda.size() < m || M.size() < m || psy.size() < m) {
    Rcpp::stop("theta_y_list slots must each have length at least m");
  }

  arma::uword p = X.n_rows;
  for (int k = 0; k < m; ++k) {
    arma::uword q_k = static_cast<arma::uword>(qVec(k));
    arma::mat y_k = Y(k);
    arma::mat lambda_k = lambda(k);
    arma::vec m_k = M(k);
    arma::mat psy_k = psy(k);

    if (y_k.n_rows != q_k || y_k.n_cols != static_cast<arma::uword>(n)) {
      Rcpp::stop("Y matrices must have q_vec[k] rows and n columns");
    }
    if (lambda_k.n_rows != p || lambda_k.n_cols != q_k) {
      Rcpp::stop("lambda matrices must have nrow(X) rows and q_vec[k] columns");
    }
    if (m_k.n_elem != p) {
      Rcpp::stop("M vectors must have length nrow(X)");
    }
    if (psy_k.n_rows != p || psy_k.n_cols != p) {
      Rcpp::stop("psy matrices must be square with dimension nrow(X)");
    }
    validate_finite_matrix(y_k, "Y");
    validate_finite_matrix(lambda_k, "lambda");
    validate_finite_matrix(psy_k, "psy");
    if (!m_k.is_finite()) {
      Rcpp::stop("M must contain only finite values");
    }
  }

  List A(m);
  arma::vec nVec(m);

  for(int k=0; k<m; ++k) {
    arma::uword q_k = static_cast<arma::uword>(qVec(k));
    arma::vec prior_diag(q_k + 1);
    prior_diag(0) = alpha1;
    prior_diag.subvec(1, q_k).fill(alpha2);
    A(k) = diagmat(prior_diag);
  };

  if (ZOneDim.n_elem != static_cast<arma::uword>(n)) {
    Rcpp::stop("length of ZOneDim must equal n");
  }
  nVec.zeros();
  arma::uvec labels(n);
  for (int i = 0; i < n; ++i) {
    double label = ZOneDim(static_cast<arma::uword>(i));
    if (!std::isfinite(label) || label < 1 || label > m || label != std::floor(label)) {
      Rcpp::stop("cluster labels must be integers in 1:m");
    }
    labels(static_cast<arma::uword>(i)) = static_cast<arma::uword>(label - 1);
  }

  List Cxxk(m);
  List Cxyk(m);
  List Cyyk(m);
  List Cytytk(m);
  List Cxtytk(m);
  List CxL1k(m);
  List Cxmyk(m);
  arma::mat sumCxmyk;
  arma::mat sumCyyk;

  for (int k=0; k<m; ++k) {
    arma::uword q_k = static_cast<arma::uword>(qVec(k));
    arma::mat Cxxkk(p, p, arma::fill::zeros);
    arma::mat Cxykk(p, q_k, arma::fill::zeros);
    arma::mat Cyykk(q_k, q_k, arma::fill::zeros);
    arma::mat Cytytkk(q_k + 1, q_k + 1, arma::fill::zeros);
    arma::mat Cxtytkk(p, q_k + 1, arma::fill::zeros);
    arma::mat CxL1kk(p, 1, arma::fill::zeros);
    arma::mat Cxmykk(p, q_k, arma::fill::zeros);
    arma::mat y_k = Y(k);
    arma::vec m_k = M(k);
    arma::mat lambda_k = lambda(k);

    for (int i = 0; i < n; ++i) {
      if (labels(static_cast<arma::uword>(i)) != static_cast<arma::uword>(k)) {
        continue;
      }
      const arma::vec x_i = X.col(i);
      const arma::vec y_i = y_k.col(i);

      nVec(k) += 1.0;
      Cxxkk += x_i * x_i.t();
      Cxykk += x_i * y_i.t();
      Cyykk += y_i * y_i.t();
      Cxtytkk.col(0) += x_i;
      Cxtytkk.cols(1, q_k) += x_i * y_i.t();
      Cytytkk(0, 0) += 1.0;
      Cytytkk.submat(0, 1, 0, q_k) += y_i.t();
      Cytytkk.submat(1, 0, q_k, 0) += y_i;
      Cytytkk.submat(1, 1, q_k, q_k) += y_i * y_i.t();
      Cxmykk += (x_i - m_k) * y_i.t();
      CxL1kk += x_i - lambda_k * y_i;
    }

    Cxxk(k) = Cxxkk;
    Cxyk(k) = Cxykk;
    Cyyk(k) = Cyykk;
    Cxtytk(k) = Cxtytkk;
    Cytytk(k) = Cytytkk;
    Cxmyk(k) = Cxmykk;
    CxL1k(k) = CxL1kk;

    arma::mat prior_cyy(q_k, q_k, arma::fill::eye);
    prior_cyy *= alpha2;
    if (k == 0) {
      sumCxmyk = Cxmykk;
      sumCyyk = Cyykk + prior_cyy;
    } else {
      sumCxmyk += Cxmykk;
      sumCyyk += Cyykk + prior_cyy;
    }
  }

  List res = Rcpp::List::create(Named("A")        = A,
                                Named("nVec")     = nVec,
                                Named("Cxxk")     = Cxxk,
                                Named("Cxyk")     = Cxyk,
                                Named("Cyyk")     = Cyyk,
                                Named("Cytytk")   = Cytytk,
                                Named("Cxtytk")   = Cxtytk,
                                Named("CxL1k")    = CxL1k,
                                Named("Cxmyk")    = Cxmyk,
                                Named("sumCxmyk") = sumCxmyk,
                                Named("sumCyyk")  = sumCyyk);

  return(res);
}
