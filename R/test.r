# library(bpgmm)
# Rcpp::sourceCpp('src/calculateCxy.cpp')
# Rcpp::sourceCpp('src/calculatePostLambdaPsy.cpp')
# # Rcpp::sourceCpp('src/dmvnorm.cpp')
# # Rcpp::sourceCpp('src/updatePostZ.cpp')
# # Rcpp::sourceCpp('src/update_Hyperparameter.cpp')
#
#
# m = 3
# n = 7
# p = 10
# q = 2
# muBar =  rnorm(p)
# dVec = c(1,1,1)
# sVec = c(1,1,1)
# qVec = c(2,2,2)
# ZOneDim = c(1,2,3,1,2,2,3)
# constraint = c(0,0,1)
# hparam <- new("Hparam", alpha1 = 3, alpha2 = 2, delta  = 3, ggamma = 4, bbeta  = 5)
# thetaYList = generatePriorThetaY(m, n, p, muBar, hparam, qVec, ZOneDim, constraint)
# X = matrix(rnorm(p * n, 0, 100), p, n)
#
#
# # a = update_Hyperparameter(m,p,q,hparam,thetaYList,dVec,sVec)
# # updateHyperparameter(m, p, q, hparam, thetaYList, dVec, sVec)
# #
# # aa = rgamma(100000, shape = 91, scale  = 0.0553164)
# #
# # mean(a)
# # var(a)
# #
# # mean(aa)
# # var(aa)
# #
# # hist(a)
# # hist(aa)
# CxyList = CalculateCxy(m, n, hparam, thetaYList, ZOneDim, qVec, X)
# a = Calculate_PostLambdaPsy(m, p, hparam, CxyList, thetaYList, qVec, constraint)
#
# alpha1 = 3
# alpha2 = 2
# delta  = 3
# ggamma = 4
# bbeta  = 5
#
# a = list()
# b = list()
# for (i in 1:1000) {
#   a[[i]] = Calculate_PostLambdaPsy(m, p, hparam, CxyList, thetaYList, qVec, constraint)
# }
#
# for (i in 1:1000) {
#   b[[i]] = R_CalculatePostLambdaPsy(alpha1, alpha2, bbeta, CxyList, thetaYList, constraint)
# }
#
# mean_lambda_cpp = 0
# sd_lambda_cpp = 0
#
# for (i in 1:1000) {
#   lambdai = a[[i]]$lambda[[1]]
#   mean_lambda_cpp = mean_lambda_cpp + mean(lambdai)
#   sd_lambda_cpp = sd_lambda_cpp + sd(lambdai)
# }
#
# mean_lambda_R = 0
# sd_lambda_R = 0
#
# for (i in 1:1000) {
#   lambdai = b[[i]]$lambda[[1]]
#   mean_lambda_R = mean_lambda_R + mean(lambdai)
#   sd_lambda_R = sd_lambda_R + sd(lambdai)
# }
#
# print(mean_lambda_cpp)
# print(mean_lambda_R)
# print(sd_lambda_cpp)
# print(sd_lambda_R)
#
#
#
# mean_psy_cpp = 0
# sd_psy_cpp = 0
#
# for (i in 1:1000) {
#   psyi = a[[i]]$psy[[1]]
#   mean_psy_cpp = mean_psy_cpp + mean(psyi)
#   sd_psy_cpp = sd_psy_cpp + sd(psyi)
# }
#
# mean_psy_R = 0
# sd_psy_R = 0
#
# for (i in 1:1000) {
#   psyi = b[[i]]$psy[[1]]
#   mean_psy_R = mean_psy_R + mean(psyi)
#   sd_psy_R = sd_psy_R + sd(psyi)
# }
#
#
# print(mean_psy_cpp)
# print(mean_psy_R)
# print(sd_psy_cpp)
# print(sd_psy_R)
#
#
#
#
#
