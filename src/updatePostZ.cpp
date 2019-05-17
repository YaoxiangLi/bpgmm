// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <iostream>
#include "utils.h"
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::S4 updatePostZ(int m,
                     int n,
                     Rcpp::S4 hparam,
                     Rcpp::S4 thetaYList,
                     arma::vec ZOneDim,
                     arma::vec qVec,
                     arma::vec constraint) {

  double alpha1 = hparam.slot("alpha1");
  double alpha2 = hparam.slot("alpha2");
  double bbeta  = hparam.slot("bbeta");

  List lambda = thetaYList.slot("lambda");
  List Y      = thetaYList.slot("Y");
  List M      = thetaYList.slot("M");
  List psy    = thetaYList.slot("psy");

  arma::vec tao;


  Rcpp::S4 res("ThetaYList");

  res.slot("tao")     = tao;
  res.slot("psy")     = psy;
  res.slot("M")       = M;
  res.slot("lambda")  = lambda;
  res.slot("Y")       = Y;

  return(res);
}

