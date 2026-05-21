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

  arma::mat Zmat = get_Z_mat(ZOneDim, m, n);
  List A(m);
  arma::vec nVec(m);

  for(int k=0; k<m; ++k) {
    nVec(k) = sum(Zmat.row(k));
    arma::vec alpha1_vec(1);

    A(k) =  diagmat(join_cols(alpha1_vec.fill(alpha1), arma::ones(qVec(k)) * alpha2));
  };

  List Cxxk;
  List Cxyk;
  List Cyyk;
  List Cytytk;
  List Cxtytk;
  List CxL1k;
  List Cxmyk;

  for (int k=0; k<m; ++k) {
    arma::mat Cxxkk;
    arma::mat Cxykk;
    arma::mat Cyykk;
    arma::mat Cytytkk;
    arma::mat Cxtytkk;
    arma::mat CxL1kk;
    arma::mat Cxmykk;
    arma::mat y_k = Y(k);
    arma::vec m_k = M(k);
    arma::mat lambda_k = lambda(k);
    arma::vec one_vec(1);

    for (int i=0; i<n; ++i) {
      if(i == 0){
        Cxxkk =  Zmat(k,i) * (X.col(i) * trans(X.col(i)));
        Cxykk =  Zmat(k,i) * (X.col(i) * trans(y_k.col(i)));
        Cyykk =  Zmat(k,i) * (y_k.col(i) * trans(y_k.col(i)));
        Cxtytkk =  Zmat(k,i) * (X.col(i) * trans(join_cols(one_vec.fill(1), y_k.col(i))));
        Cytytkk =  Zmat(k,i) * (join_cols(one_vec.fill(1), y_k.col(i)) * trans(join_cols(one_vec.fill(1), y_k.col(i))));
        Cxmykk =  Zmat(k,i) * ((X.col(i) - m_k) * trans(y_k.col(i)));
        CxL1kk = Zmat(k,i) * (X.col(i) -  (lambda_k * y_k.col(i)));

      }else{
        Cxxkk = Cxxkk + Zmat(k,i) * (X.col(i) * trans(X.col(i)));
        Cxykk = Cxykk + Zmat(k,i) * (X.col(i) * trans(y_k.col(i)));
        Cyykk = Cyykk +  Zmat(k,i) * (y_k.col(i) * trans(y_k.col(i)));
        Cxtytkk = Cxtytkk + Zmat(k,i) * (X.col(i) * trans(join_cols(one_vec.fill(1), y_k.col(i))));
        Cytytkk = Cytytkk + Zmat(k,i) * (join_cols(one_vec.fill(1), y_k.col(i)) * trans(join_cols(one_vec.fill(1), y_k.col(i))));
        Cxmykk = Cxmykk + Zmat(k,i) * ((X.col(i) - m_k) * trans(y_k.col(i)));
        CxL1kk = CxL1kk + Zmat(k,i) * (X.col(i) -  (lambda_k * y_k.col(i)));
      }
    }
    Cxxk.push_back(Cxxkk);
    Cxyk.push_back(Cxykk);
    Cyyk.push_back(Cyykk);
    Cxtytk.push_back(Cxtytkk);
    Cytytk.push_back(Cytytkk);
    Cxmyk.push_back(Cxmykk);
    CxL1k.push_back(CxL1kk);
  }

  arma::mat sumCxmyk;
  arma::mat sumCyyk;

  for (int k=0; k<m; ++k) {
      arma::mat cxmyk_k = Cxmyk[k];
      arma::mat cyyk_k = Cyyk[k];
    if (k == 0) {
      sumCxmyk = cxmyk_k;
      sumCyyk  = cyyk_k + diagmat( arma::ones(qVec(k))) * alpha2;
    } else {
      sumCxmyk = sumCxmyk + cxmyk_k;
      sumCyyk  = sumCyyk + cyyk_k + diagmat(arma::ones(qVec(k))) * alpha2;
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
