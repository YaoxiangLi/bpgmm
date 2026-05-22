// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <cmath>
#include "utils.h"

namespace {

bool is_active(const arma::vec &clus_ind, int k) {
  const double value = clus_ind(static_cast<arma::uword>(k));
  if (!std::isfinite(value) || (value != 0.0 && value != 1.0)) {
    Rcpp::stop("clusInd entries must be 0/1");
  }
  return value == 1.0;
}

arma::mat list_matrix(const Rcpp::List &values,
                      int index,
                      arma::uword expected_rows,
                      arma::uword expected_cols,
                      const char* name) {
  arma::mat out = Rcpp::as<arma::mat>(values[index]);
  if (out.n_rows != expected_rows || out.n_cols != expected_cols || !out.is_finite()) {
    Rcpp::stop("%s entries have incompatible dimensions or non-finite values", name);
  }
  return out;
}

arma::mat list_matrix_with_rows(const Rcpp::List &values,
                                int index,
                                arma::uword expected_rows,
                                const char* name) {
  arma::mat out = Rcpp::as<arma::mat>(values[index]);
  if (out.n_rows != expected_rows || out.n_cols < 1 || !out.is_finite()) {
    Rcpp::stop("%s entries have incompatible dimensions or non-finite values", name);
  }
  return out;
}

double log_gamma_rate(double x, double shape, double rate) {
  if (!std::isfinite(x) || x <= 0.0) {
    Rcpp::stop("gamma density inputs must be positive finite values");
  }
  return R::dgamma(x, shape, 1.0 / rate, true);
}

double log_mvn_column(const arma::vec &x,
                      const arma::mat &sigma) {
  arma::mat row_x = x.t();
  arma::rowvec mean = arma::zeros<arma::rowvec>(x.n_elem);
  return dmvnrm_arma(row_x, mean, sigma, true)(0);
}

arma::mat harmonic_psi_average(const Rcpp::List &psy, int p, int m, const arma::vec &clus_ind) {
  arma::mat precision_sum = arma::zeros<arma::mat>(p, p);
  for (int k = 0; k < m; ++k) {
    if (is_active(clus_ind, k)) {
      precision_sum += arma::inv(list_matrix(psy, k, p, p, "psy"));
    }
  }
  return arma::inv(precision_sum / static_cast<double>(m));
}

} // namespace

// [[Rcpp::export]]
double Evaluate_PriorPsi(Rcpp::List psy,
                         int p,
                         int m,
                         double delta,
                         double bbeta,
                         arma::vec constraint,
                         arma::vec clusInd) {
  validate_positive_int(p, "p");
  validate_positive_int(m, "m");
  validate_constraint_vec(constraint);
  if (!std::isfinite(delta) || delta <= 0.0 || !std::isfinite(bbeta) || bbeta <= 0.0) {
    Rcpp::stop("delta and bbeta must be positive finite values");
  }
  if (psy.size() < m || clusInd.n_elem < static_cast<arma::uword>(m)) {
    Rcpp::stop("psy and clusInd must contain at least m entries");
  }

  double value = 0.0;
  const bool common_psi = constraint(1) == 1.0;
  const bool isotropic = constraint(2) == 1.0;

  for (int k = 0; k < m; ++k) {
    if (!is_active(clusInd, k) || (common_psi && k != 0)) {
      continue;
    }

    arma::mat psi_k = list_matrix(psy, k, p, p, "psy");
    if (isotropic) {
      value += log_gamma_rate(1.0 / psi_k(0, 0), delta, bbeta);
    } else {
      arma::vec inv_diag = 1.0 / psi_k.diag();
      for (arma::uword j = 0; j < inv_diag.n_elem; ++j) {
        value += log_gamma_rate(inv_diag(j), delta, bbeta);
      }
    }
  }

  return value;
}

// [[Rcpp::export]]
double Evaluate_PriorLambda(int p,
                            int m,
                            double alpha2,
                            arma::vec qVec,
                            Rcpp::List psy,
                            Rcpp::List lambda,
                            arma::vec constraint,
                            arma::vec clusInd) {
  validate_positive_int(p, "p");
  validate_positive_int(m, "m");
  validate_q_vec(qVec, m);
  validate_constraint_vec(constraint);
  if (!std::isfinite(alpha2) || alpha2 <= 0.0) {
    Rcpp::stop("alpha2 must be a positive finite value");
  }
  if (psy.size() < m || lambda.size() < m || clusInd.n_elem < static_cast<arma::uword>(m)) {
    Rcpp::stop("psy, lambda, and clusInd must contain at least m entries");
  }

  const bool common_lambda = constraint(0) == 1.0;
  const bool common_psi = constraint(1) == 1.0;
  arma::mat shared_psi_average;
  if (common_lambda && !common_psi) {
    shared_psi_average = harmonic_psi_average(psy, p, m, clusInd);
  }

  double value = 0.0;
  for (int k = 0; k < m; ++k) {
    if (!is_active(clusInd, k) || (common_lambda && k != 0)) {
      continue;
    }

    arma::mat lambda_k = list_matrix_with_rows(lambda, k, p, "lambda");
    arma::mat sigma = (1.0 / alpha2) *
      ((common_lambda && !common_psi) ? shared_psi_average : list_matrix(psy, k, p, p, "psy"));

    for (arma::uword j = 0; j < lambda_k.n_cols; ++j) {
      value += log_mvn_column(lambda_k.col(j), sigma);
    }
  }

  return value;
}
