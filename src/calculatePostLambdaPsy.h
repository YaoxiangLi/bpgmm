#ifndef __BPGMM_CALCULATEPOSTLAMBDAPSY__
#define __BPGMM_CALCULATEPOSTLAMBDAPSY__

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
