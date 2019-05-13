# library(bpgmm)
# Rcpp::sourceCpp('src/calculatePostLambdaPsy.cpp')
# m = 3
# n = 7
# p = 10
# muBar =  rnorm(p)
# qVec = c(2,2,2)
# ZOneDim = c(1,2,3,1,2,2,3)
# constraint = c(1,1,1)
# hparam <- new("Hparam",alpha1 = 3,
#               alpha2 = 2,
#               delta  = 3,
#               ggamma = 4,
#               bbeta  = 5)
# thetaYList = generatePriorThetaY(m, n, p, muBar, hparam, qVec, ZOneDim, constraint)
# X = matrix(rnorm(p * n, 0, 1), p, n)
#
# cxyList = CalculateCxy(m, n, hparam, thetaYList, ZOneDim, qVec, X)
# r_cxyList = r_CalculateCxy(m, n, hparam, thetaYList, ZOneDim, qVec, X)
#
# CalculatePostLambdaPsy(hparam, cxyList, thetaYList, constraint)
