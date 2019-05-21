// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <iostream>
#include "utils.h"
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::NumericVector updatePostZ(
                     arma::mat X,
                     int m,
                     int n,
                     Rcpp::S4 thetaYList ){

  // RInside R(1, "2");
  // Obtaining namespace of rmvnorm function
  Rcpp::Environment mvtnorm = Rcpp::Environment::namespace_env("mvtnorm");
  // Picking up rmvnorm() function from rmvnorm package
  Rcpp::Function dmvnorm = mvtnorm["dmvnorm"];

  Rcpp::NumericVector ZOneDim;

  List lambda = thetaYList.slot("lambda");
  List Y      = thetaYList.slot("Y");
  List M      = thetaYList.slot("M");
  List psy    = thetaYList.slot("psy");


  arma::mat pMat(m,n);
  arma::mat dMat(m,n);
//   for(k in 1:m){
//     for(i in 1:n){
//       dMat[k,i] = dmvnorm(X[,i], mean = M[[k]], sigma = psy[[k]] +  lambda[[k]]%*%t(lambda[[k]])
//                             ,log = T)
//     }
//   }

  // double a[] = {1,2};
  // Rcpp::NumericVector aa = Rcpp::as<Rcpp::NumericVector> (wrap(a));
  // Rcpp::NumericVector temp= dmvnorm(Named("x", Rcpp::as<Rcpp::NumericVector> (wrap(a))),
  //                                 Named("mean", Rcpp::as<Rcpp::NumericVector> (wrap(a))));
  // std::cout << temp << std::endl;




  for(int k = 0; k < m; k++ ){
    for(int i = 0; i < n; i++){

      arma::mat Xi = X.col(i);
      arma::vec Mk = M(k);
      arma::mat lambdak = lambda(k);
      arma::mat psyk = psy(k);
      arma::mat var = psyk +  lambdak * trans(lambdak);


      // NumericVector Mk = M(k);
      // arma::vec Mk2 = Rcpp::as<arma::vec> (wrap(Mk));
      // Rcpp::NumericVector Mk3 = Rcpp::as<Rcpp::NumericVector> (wrap(Mk2));
      //
      // Rcpp::NumericVector temp= dmvnorm(Named("x", Mk),Named("mean", Mk));
      // Rcpp::NumericVector temp2= dmvnorm(Named("x", Mk2),Named("mean", Mk2));
      // Rcpp::NumericVector temp3= dmvnorm(Named("x", Mk3),Named("mean", Mk3));


      arma::vec temp =  dmvnrm_arma(Xi, Mk, var,true);
       std::cout << "temp" << std::endl << temp << std::endl;
      // std::cout << "MK" << std::endl << Mk << std::endl;
      // std::cout << "MK1" << std::endl <<Mk2 << std::endl;
      // std::cout << "MK2" << std::endl << Mk3<< std::endl;

      // Rcpp::NumericVector temp= dmvnorm(Named("x", Mk3),
      //                                    Named("mean", Mk3));
      // std::cout << temp << std::endl;
      // Rcpp::NumericVector temp2= dmvnorm(Named("x", Mk),
      //                                    Named("mean", Mk));


      // std::cout << "mean" << std::endl << Rcpp::as<Rcpp::NumericVector>(wrap(Mk)) << std::endl;
      // std::cout << "var" << std::endl << Rcpp::as<Rcpp::NumericMatrix>(wrap(var)) << std::endl;
      // std::cout << "Xi" << std::endl << Rcpp::as<Rcpp::NumericVector>(wrap(Xi)) << std::endl;
      // dMat(k,i)
      // Rcpp::NumericVector mean1 = Rcpp::as<Rcpp::NumericVector>(wrap(Mk));
      // Rcpp::NumericVector Xi1 = Rcpp::as<Rcpp::NumericVector>(wrap(Xi));
      // Rcpp::NumericVector temp= dmvnorm(Named("x", Xi1),Named("mean", mean1));

      // std::cout << "mean" << std::endl << Mk << std::endl;
      // std::cout << "Xi" << std::endl << Xi << std::endl;

      // Rcpp::NumericVector temp2= dmvnorm(Named("x", Rcpp::as<arma::vec>(wrap(Xi))),
      //                                    Named("mean", Rcpp::as<arma::vec>(wrap(Mk))));

      // Named("sigma",Rcpp::as<Rcpp::NumericMatrix>(wrap(var))));
      // std::cout << temp2 << std::endl;
      // Rcpp::NumericVector dki= dmvnorm(Named("x", Xi),
      //                   Named("mean", Mk),
      //                   Named("log", "T"));
      // Rcpp::NumericVector temp = Rcpp::as(dki);
      // dMat(k,i)
      // arma::vec temp  = Rcpp::as(dki);
      // std::cout << dki << std::endl;
    }
  }
 // std::cout << dMat << std::endl;


//   for(k in 1:m){
//     dMat[k,] = dMat[k,] + log(tao[k])
//   }
//
//   for(i in 1:n){
//     for(k in 1:m){
//       pMat[k,i] = calculateRatio(dMat[k,i], dMat[,i])
//     }
//   }
//
//   ZOneDim = c()
//     for(i in 1:n){
//       tempProb = as.numeric(pMat[,i])
//       ZOneDim[i] = sample(x = 1:m, size = 1, prob = tempProb)
//     }
//     ZOneDim


  arma::vec tao;




  // Rcpp::S4 res("ThetaYList");
  // res.slot("tao")     = tao;
  // res.slot("psy")     = psy;
  // res.slot("M")       = M;
  // res.slot("lambda")  = lambda;
  // res.slot("Y")       = Y;

  return(ZOneDim);
}

