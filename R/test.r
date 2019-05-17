library(bpgmm)
Rcpp::sourceCpp('src/calculateCxy.cpp')
Rcpp::sourceCpp('src/calculatePostLambdaPsy.cpp')


m = 3
n = 7
p = 10
muBar =  rnorm(p)
qVec = c(2,2,2)
ZOneDim = c(1,2,3,1,2,2,3)
constraint = c(1,1,0)
hparam <- new("Hparam", alpha1 = 3, alpha2 = 2, delta  = 3, ggamma = 4, bbeta  = 5)
thetaYList = generatePriorThetaY(m, n, p, muBar, hparam, qVec, ZOneDim, constraint)
X = matrix(rnorm(p * n, 0, 1), p, n)

CxyList = CalculateCxy(m, n, hparam, thetaYList, ZOneDim, qVec, X)
CalculatePostLambdaPsy(m, p, hparam, CxyList, thetaYList, qVec, constraint)

k = 1
Cxxk = CxyList$Cxxk
Cxtytk = CxyList$Cxtytk
thetaYList@psy
matrix(c)


diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
     + tildaLambda[[k]]%*%(Cytytk[[k]] + A[[k]])%*%t(tildaLambda[[k]] )
     + 2 * bbeta * diag(rep(1,p)))




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

