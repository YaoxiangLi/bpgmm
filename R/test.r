Rcpp::sourceCpp('src/calculateCxy.cpp')
m = 3
n = 7
p = 10
muBar =  rnorm(p)
qVec = c(2,2,2)
ZOneDim = c(1,2,3,1,2,2,3)
constraint = c(1,1,1)
hparam <- new("Hparam",alpha1 = 3,
              alpha2 = 2,
              delta  = 3,
              ggamma = 4,
              bbeta  = 5)
thetaYList = generatePriorThetaY(m, n, p, muBar, hparam, qVec, ZOneDim, constraint)
X = matrix(rnorm(p * n, 0, 1), p, n)
res = CalculateCxy(m, n, hparam, thetaYList, ZOneDim, qVec, X)

r_res = r_CalculateCxy(m, n, hparam, thetaYList, ZOneDim, qVec, X)

microbenchmarkCore::microbenchmark(CalculateCxy(m, n, hparam, thetaYList, ZOneDim, qVec, X))
microbenchmarkCore::microbenchmark(r_CalculateCxy(m, n, hparam, thetaYList, ZOneDim, qVec, X))

testthat::expect_equal(r_res$A, res$A)
testthat::expect_equal(r_res$nVec, res$nVec)
testthat::expect_equal(r_res$Cxxk, res$Cxxk)
testthat::expect_equal(r_res$Cxyk, res$Cxyk)
testthat::expect_equal(r_res$Cyyk, res$Cyyk)
testthat::expect_equal(r_res$Cytytk, res$Cytytk)
testthat::expect_equal(r_res$Cxtytk[1], res$Cxtytk[1])
testthat::expect_equal(r_res$CxL1k, res$CxL1k)
testthat::expect_equal(r_res$Cxmyk, res$Cxmyk)
testthat::expect_equal(r_res$sumCxmyk, res$sumCxmyk)
testthat::expect_equal(r_res$sumCyyk, res$sumCyyk)



r_res$Cytytk[[1]]
res$Cytytk[[1]]
r_res$Cxtytk
res$Cxtytk



