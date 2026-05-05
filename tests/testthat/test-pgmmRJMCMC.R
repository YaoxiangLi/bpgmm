test_that("pgmmRJMCMC default constraint is a numeric vector", {
  x <- matrix(rnorm(12), nrow = 2)

  expect_error(
    pgmmRJMCMC(x, mInit = 1, mVec = c(1, 2), qnew = 1, burn = 0, niter = 1),
    NA
  )
})
