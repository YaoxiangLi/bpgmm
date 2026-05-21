#ifndef __BPGMM_UTILS__
#define __BPGMM_UTILS__

#include <RcppArmadillo.h>

arma::vec dmvnrm_arma(const arma::mat& x,
                      const arma::rowvec& mean,
                      const arma::mat& sigma,
                      bool logd);


arma::mat get_Z_mat(const arma::vec& ZOneDim, int m, int n);

double calculate_Ratio(double logDeno, const arma::vec& logNume);

void validate_positive_int(int value, const char* name);

void validate_finite_matrix(const arma::mat& x, const char* name);

void validate_q_vec(const arma::vec& q_vec, int m);

void validate_constraint_vec(const arma::vec& constraint);

void validate_positive_finite_vec(const arma::vec& x,
                                  arma::uword expected_length,
                                  const char* name);

double get_positive_finite_slot(Rcpp::S4 obj, const char* slot_name);

Rcpp::IntegerVector update_PostZ(arma::mat X,
                  int m,
                  int n,
                  Rcpp::S4 thetaYList,
                  arma::vec tao);

Rcpp::List Calculate_Cxy(int m, int n, Rcpp::S4 hparam, Rcpp::S4 thetaYList,
                  arma::vec ZOneDim,
                  arma::vec qVec,
                  arma::mat X);

Rcpp::List Calculate_PostLambdaPsy(int m,
                                   int p,
                                   Rcpp::S4 hparam,
                                   Rcpp::List CxyList,
                                   Rcpp::S4 thetaYList,
                                   arma::vec qVec,
                                   arma::vec constraint);

Rcpp::S4 update_Hyperparameter(int m,
                               int p,
                               int q,
                               Rcpp::S4 hparam,
                               Rcpp::S4 thetaYList,
                               arma::vec dVec,
                               arma::vec sVec);


#endif
