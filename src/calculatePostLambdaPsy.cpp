// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <iostream>
#include "utils.h"
using namespace Rcpp;



// [[Rcpp::export]]
Rcpp::List CalculatePostLambdaPsy(Rcpp::S4 hparam,
                                  Rcpp::List CxyList,
                                  Rcpp::S4 thetaYList,
                                  arma::vec constraint) {

  double alpha1 = hparam.slot("alpha1");
  double alpha2 = hparam.slot("alpha2");
  double bbeta  = hparam.slot("bbeta");

  List Cxxk      = CxyList["Cxxk"];
  List Cxyk      = CxyList["Cxyk"];
  List Cyyk      = CxyList["Cyyk"];
  List Cytytk    = CxyList["Cytytk"];
  List Cxtytk    = CxyList["Cxtytk"];
  List CxL1k     = CxyList["CxL1k"];
  List Cxmyk     = CxyList["Cxmyk"];

  arma::mat sumCxmyk  = CxyList["sumCxmyk"];
  arma::mat sumCyyk   = CxyList["sumCyyk"];

  List A         = CxyList["A"];
  arma::vec nVec = CxyList["nVec"];

  List M      = thetaYList.slot("M");
  List psy    = thetaYList.slot("psy");

  List lambda;

  std::cout << constraint << std::endl;



  List res = Rcpp::List::create(Named("lambda") = lambda,
                                Named("psy")    = psy);

  return(res);
}




