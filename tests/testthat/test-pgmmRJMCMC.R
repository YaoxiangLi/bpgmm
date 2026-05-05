test_that("pgmm_rjmcmc default constraint is a numeric vector", {
  x <- matrix(rnorm(12), nrow = 2)

  expect_error(
    pgmm_rjmcmc(x, mInit = 1, mVec = c(1, 2), qnew = 1, burn = 0, niter = 1),
    NA
  )
})

test_that("pgmmRJMCMC is a deprecated compatibility wrapper", {
  x <- matrix(rnorm(12), nrow = 2)

  expect_warning(
    expect_error(
      pgmmRJMCMC(x, mInit = 1, mVec = c(1, 2), qnew = 1, burn = 0, niter = 1),
      NA
    ),
    "deprecated"
  )
})
