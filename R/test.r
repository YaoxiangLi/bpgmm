library(bpgmm)
Rcpp::sourceCpp('src/calculateCxy.cpp')
Rcpp::sourceCpp('src/calculatePostLambdaPsy.cpp')
# Rcpp::sourceCpp('src/dmvnorm.cpp')
# Rcpp::sourceCpp('src/updatePostZ.cpp')
m = 3
n = 7
p = 10
muBar =  rnorm(p)
qVec = c(2,2,2)
ZOneDim = c(1,2,3,1,2,2,3)
constraint = c(0,0,1)
hparam <- new("Hparam", alpha1 = 3, alpha2 = 2, delta  = 3, ggamma = 4, bbeta  = 5)
thetaYList = generatePriorThetaY(m, n, p, muBar, hparam, qVec, ZOneDim, constraint)
X = matrix(rnorm(p * n, 0, 1), p, n)

# updatePostZ(m, n, thetaYList)


CxyList = CalculateCxy(m, n, hparam, thetaYList, ZOneDim, qVec, X)
CalculatePostLambdaPsy(m, p, hparam, CxyList, thetaYList, qVec, constraint)




k = 1
Cxxk = CxyList$Cxxk
Cxmyk = CxyList$Cxmyk
Cyyk = CxyList$Cyyk

Cxtytk = CxyList$Cxtytk
thetaYList@psy
matrix(c)






microbenchmarkCore::microbenchmark(CalculatePostLambdaPsy(m, hparam, cxyList, thetaYList, constraint))
constraint = c(1,1,0)
microbenchmarkCore::microbenchmark(CalculatePostLambdaPsy(m, hparam, cxyList, thetaYList, constraint))
constraint = c(1,0,1)
microbenchmarkCore::microbenchmark(CalculatePostLambdaPsy(m, hparam, cxyList, thetaYList, constraint))
constraint = c(1,0,0)
microbenchmarkCore::microbenchmark(CalculatePostLambdaPsy(m, hparam, cxyList, thetaYList, constraint))
constraint = c(0,1,1)
microbenchmarkCore::microbenchmark(CalculatePostLambdaPsy(m, hparam, cxyList, thetaYList, constraint))
constraint = c(0,1,0)
microbenchmarkCore::microbenchmark(CalculatePostLambdaPsy(m, hparam, cxyList, thetaYList, constraint))
constraint = c(0,0,1)
microbenchmarkCore::microbenchmark(CalculatePostLambdaPsy(m, hparam, cxyList, thetaYList, constraint))
constraint = c(0,0,0)
microbenchmarkCore::microbenchmark(CalculatePostLambdaPsy(m, hparam, cxyList, thetaYList, constraint))

