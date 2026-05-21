#ifndef BPGMM_CALCULATE_POST_LAMBDA_PSY_H
#define BPGMM_CALCULATE_POST_LAMBDA_PSY_H

#include <RcppArmadillo.h>
#include "utils.h"

Rcpp::List Calculate_PostLambdaPsy(int m,
                                   int p,
                                   Rcpp::S4 hparam,
                                   Rcpp::List CxyList,
                                   Rcpp::S4 thetaYList,
                                   arma::vec qVec,
                                   arma::vec constraint);


#endif
