// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::plugins(cpp11)]]
// #include <RcppArmadillo.h>
// #include <iostream>
// #include "utils.h"
// using namespace Rcpp;
//
//
// // [[Rcpp::export]]
// Rcpp::NumericVector updatePostZ(
//                      int m,
//                      int n,
//                      Rcpp::S4 thetaYList ){
//
//   // Obtaining namespace of rmvnorm function
//   Rcpp::Environment mvtnorm = Rcpp::Environment::namespace_env("mvtnorm");
//   // Picking up rmvnorm() function from rmvnorm package
//   Rcpp::Function rmvnorm = mvtnorm["rmvnorm"];
//
//   Rcpp::NumericVector ZOneDim;
//
//   List lambda = thetaYList.slot("lambda");
//   List Y      = thetaYList.slot("Y");
//   List M      = thetaYList.slot("M");
//   List psy    = thetaYList.slot("psy");
//
//
//   Rcpp::NumericMatrix pMat(m,n);
//   Rcpp::NumericMatrix dMat(m,n);
// //   for(k in 1:m){
// //     for(i in 1:n){
// //       dMat[k,i] = dmvnorm(X[,i], mean = M[[k]], sigma = psy[[k]] +  lambda[[k]]%*%t(lambda[[k]])
// //                             ,log = T)
// //     }
// //   }
//
//   for(int k = 0; k < m; k++ ){
//     for(int i = 0; i < n; i++){
//
//       Rcpp::NumericVector Mk = M[k];
//       Rcpp::NumericMatrix lambdak = lambda[k] ;
//
//       Rcpp::NumericMatrix lambdakk = lambdak * transpose(lambdak);
//
//       std::cout << "mean" << std::endl << Mk << std::endl;;
//       std::cout << "var" << std::endl << lambdakk << std::endl;;
//       dMat(k,i) = k+i;
//
//     }
//   }
// //  std::cout << dMat << std::endl;
//
//
// //   for(k in 1:m){
// //     dMat[k,] = dMat[k,] + log(tao[k])
// //   }
// //
// //   for(i in 1:n){
// //     for(k in 1:m){
// //       pMat[k,i] = calculateRatio(dMat[k,i], dMat[,i])
// //     }
// //   }
// //
// //   ZOneDim = c()
// //     for(i in 1:n){
// //       tempProb = as.numeric(pMat[,i])
// //       ZOneDim[i] = sample(x = 1:m, size = 1, prob = tempProb)
// //     }
// //     ZOneDim
//
//
//   arma::vec tao;
//
//
//
//
//   // Rcpp::S4 res("ThetaYList");
//   // res.slot("tao")     = tao;
//   // res.slot("psy")     = psy;
//   // res.slot("M")       = M;
//   // res.slot("lambda")  = lambda;
//   // res.slot("Y")       = Y;
//
//   return(ZOneDim);
// }
//
