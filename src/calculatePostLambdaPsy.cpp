// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <iostream>
#include "utils.h"
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
Rcpp::List CalculatePostLambdaPsy(int m,
                                  int p,
                                  Rcpp::S4 hparam,
                                  Rcpp::List CxyList,
                                  Rcpp::S4 thetaYList,
                                  arma::vec qVec,
                                  arma::vec constraint) {

  // double alpha1 = hparam.slot("alpha1");
  // double alpha2 = hparam.slot("alpha2");
  double bbeta  = hparam.slot("bbeta");
  double delta  = hparam.slot("delta");

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

  List lambda(m);

  // std::cout << constraint << std::endl;
  // std::cout << constraint[0] << std::endl;
  // std::cout << constraint[1] << std::endl;
  // std::cout << constraint[2] << std::endl;


  // Obtaining namespace of rmvnorm function
  Rcpp::Environment mvtnorm = Rcpp::Environment::namespace_env("mvtnorm");
  // Picking up rmvnorm() function from rmvnorm package
  Rcpp::Function rmvnorm = mvtnorm["rmvnorm"];


  Rcpp::Environment base = Rcpp::Environment::namespace_env("base");
  Rcpp::Environment stats = Rcpp::Environment::namespace_env("stats");
  Rcpp::Function kronecker = base["kronecker"];
  Rcpp::Function c = base["c"];
  Rcpp::Function matrix = base["matrix"];
  Rcpp::Function t = base["t"];
  Rcpp::Function diag = base["diag"];
  Rcpp::Function rgamma = stats["rgamma"];


  // Rcpp::NumericVector test_norm = c(sumCxmyk * sumCyyk.i());
  // test_norm = c(sumCxmyk * sumCyyk.i());
  // std::cout  << test_norm << std::endl;


  // Rcpp::NumericVector test_norm1 = kronecker(sumCyyk.i(), psy[1]);
  // std::cout  << test_norm << std::endl;

  // std::cout << sumCyyk << sumCyyk.i() << std::endl;
  // arma::vec test_norm;
  // test_norm = rmvnorm(_["n"] = 1, _["mean"]=arma::vec(1), _["sigma"]=arma::vec(1));

  if (constraint[0] == 1 & constraint[1] == 1 & constraint[2] == 1) {
    std::cout << "Model 1" << std::endl;
    // Model 1

    for (int k=0; k<m; ++k) {
      if (k == 0) {

        Rcpp::NumericVector mean_vec = c(sumCxmyk * sumCyyk.i());
        // std::cout << sumCxmyk << std::endl;
        // std::cout << sumCyyk.i() << std::endl;
        // std::cout << mean_vec << std::endl;

        Rcpp::NumericMatrix sigma_mat = kronecker(sumCyyk.i(), psy[k]);
        // std::cout << sigma_mat << std::endl;


        lambda[k] = rmvnorm(Named("n", 1),
                            Named("mean", mean_vec),
                            Named("sigma", sigma_mat));
        Rcpp::NumericMatrix lambdak = lambda[k];
        // std::cout << lambdak << std::endl;
        lambda[k] = matrix(lambdak, p, qVec[k]);

      } else{
        lambda[k] = lambda[0];
      }
    }

    // post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k

    List tildaLambda(m);
    for (int k=0; k<m; ++k) {
      Rcpp::NumericMatrix lambda_k = lambda[k];
      Rcpp::NumericVector m_k = M[k];
      // std::cout << "lambda_k" << lambda_k << std::endl;
      // std::cout << "m_k: " << m_k << std::endl;

      arma::vec m_ka = m_k;
      // std::cout << "m_ka: "  << std::endl << m_ka << std::endl;
      arma::mat lambda_ka = Rcpp::as<arma::mat>(lambda_k);
      // std::cout << "lambda_ka: "  << std::endl << lambda_ka << std::endl;
      lambda_ka.insert_cols(0, m_ka);
      // std::cout << "lambda_ka: " << std::endl << lambda_ka << std::endl;
      tildaLambda[k] = lambda_ka;
    }

    // Post psy, psy was defined in model before, here for clear?

    List post_psy(m);

    double shapePara = 0;
    arma::vec ratePara_vec(p, arma::fill::zeros);
    // std::cout << "ratePara_vec: "  << std::endl << ratePara_vec << std::endl;

    for (int k=0; k<m; ++k) {

      // double nVec_k = nVec[k];
      // double qVec_k = qVec[k];


      // std::cout << "nVec[k]: "  << std::endl << nVec[k] << std::endl;
      // std::cout << "qVec[k]: "  << std::endl << qVec[k] << std::endl;
      // std::cout << "delta: "  << std::endl << delta << std::endl;


      shapePara += p/2 * (nVec[k] + qVec[k] + 2 * delta - 1);
      // std::cout << "shapePara: "  << std::endl << shapePara << std::endl;
      // ratePara_vec += 1；

      Rcpp::NumericMatrix Cxxk_k = Cxxk[k];
      Rcpp::NumericMatrix Cxtytk_k = Cxtytk[k];
      Rcpp::NumericMatrix Cytytk_k = Cytytk[k];

      Rcpp::NumericMatrix tildaLambda_k = tildaLambda[k];
      Rcpp::NumericMatrix A_k = A[k];

      // std::cout << "Cxxk[k]: " << std::endl << Cxxk_k << std::endl;
      // std::cout << "Cxtytk_k: " << std::endl << Cxtytk_k << std::endl;
      // std::cout << "tildaLambda_k: " << std::endl << tildaLambda_k << std::endl;

      arma::mat Cxxk_ka = Rcpp::as<arma::mat>(Cxxk_k);
      arma::mat Cxtytk_ka = Rcpp::as<arma::mat>(Cxtytk_k);
      arma::mat Cytytk_ka = Rcpp::as<arma::mat>(Cytytk_k);

      arma::mat tildaLambda_ka = Rcpp::as<arma::mat>(tildaLambda_k);
      arma::mat A_ka = Rcpp::as<arma::mat>(A_k);
      arma::mat bbeta_eye(p, p, arma::fill::eye);

      bbeta_eye = 2 * bbeta * bbeta_eye;
      // std::cout << "Cxxk[k]: (p * p)" << std::endl << Cxxk_k << std::endl;
      // std::cout << "Cxtytk_ka: (p * m)" << std::endl << Cxtytk_ka << std::endl;
      // std::cout << "tildaLambda_ka: (m * p)" << std::endl << trans(tildaLambda_ka) << std::endl;
      // std::cout << "A_ka" << std::endl <<  A_ka << std::endl;
      // std::cout << "bbeta_eye" << std::endl <<  bbeta_eye << std::endl;

      arma::mat ratePara_k = Cxxk_ka - 2 * Cxtytk_ka * trans(tildaLambda_ka) + tildaLambda_ka * (Cytytk_ka + A_ka) * trans(tildaLambda_ka) + bbeta_eye;
      // ratePara_k = 1/2 * arma::diagvec(ratePara_k);
      ratePara_vec += arma::diagvec(ratePara_k) / 2;
      // std::cout << "ratePara_k: " << std::endl << ratePara_k << std::endl;

      // std::cout << "ratePara_k: " << std::endl << arma::diagvec(ratePara_k) / 2 << std::endl;
      // std::cout << "ratePara_vec: " << std::endl << ratePara_vec << std::endl;

      // std::cout << "test: (p * p)" << std::endl << arma::diagvec(test) << std::endl;


    }
    shapePara += 1;
    double ratePara = arma::sum(ratePara_vec);
    double scalePara = 1 / ratePara;
    // std::cout << "shapePara: " << std::endl << shapePara << std::endl;
    // std::cout << "ratePara: " << std::endl << ratePara << std::endl;
    // std::cout << "scalePara: " << std::endl << scalePara << std::endl;
    double invpsy = sum(arma::randg( 1, distr_param(shapePara, scalePara) ));
    // std::cout << "invpsy: " << std::endl << invpsy << std::endl;

    // double invpsy = rgamma(Named("n", 1), Named("shape", shapePara), Named("rate", ratePara));
    // std::cout << "invpsy: " << std::endl << invpsy << std::endl;

    for (int k=0; k<m; ++k) {
      arma::mat post_psy_eye(p, p, arma::fill::eye);
      // std::cout << "post_psy_eye: " << std::endl << post_psy_eye << std::endl;
      // std::cout << "post_psy_eye: " << std::endl << 1/invpsy * post_psy_eye << std::endl;

      post_psy(k) = 1/invpsy * post_psy_eye;
    }
    List res = Rcpp::List::create(Named("lambda") = lambda,
                                  Named("psy")    = post_psy);

    return(res);

  } else if (constraint[0] == 1 & constraint[1] == 1 & constraint[2] == 0) {
    std::cout << "Model 2" << std::endl;
  } else if (constraint[0] == 1 & constraint[1] == 0 & constraint[2] == 1) {
    std::cout << "Model 3" << std::endl;
  } else if (constraint[0] == 1 & constraint[1] == 0 & constraint[2] == 0) {
    std::cout << "Model 4" << std::endl;
  } else if (constraint[0] == 0 & constraint[1] == 1 & constraint[2] == 1) {
    std::cout << "Model 5" << std::endl;
  } else if (constraint[0] == 0 & constraint[1] == 1 & constraint[2] == 0) {
    std::cout << "Model 6" << std::endl;
  } else if (constraint[0] == 0 & constraint[1] == 0 & constraint[2] == 1) {
    std::cout << "Model 7" << std::endl;
  } else if (constraint[0] == 0 & constraint[1] == 0 & constraint[2] == 0) {
    std::cout << "Model 8" << std::endl;
  }


  // Test Results


  // True Results
  // List res = Rcpp::List::create(Named("lambda") = lambda,
  //                               Named("psy")    = psy);

  // return(res);
}




