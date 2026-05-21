test_that("pgmm_rjmcmc default constraint is a numeric vector", {
  x <- matrix(rnorm(12), nrow = 2)

  expect_error(
    pgmm_rjmcmc(x, mInit = 1, mVec = c(1, 2), qnew = 1, burn = 0, niter = 1, verbose = FALSE),
    NA
  )
})

test_that("pgmm_rjmcmc accepts paper model labels for constraints", {
  x <- matrix(rnorm(12), nrow = 2)

  fit <- pgmm_rjmcmc(
    x,
    mInit = 1,
    mVec = c(1, 2),
    qnew = 1,
    burn = 0,
    niter = 1,
    constraint = "UUU",
    verbose = FALSE
  )

  expect_equal(fit$constraintList[[1]], model_to_constraint("UUU"))
})

test_that("pgmmRJMCMC is a deprecated compatibility wrapper", {
  x <- matrix(rnorm(12), nrow = 2)

  expect_warning(
    expect_error(
      pgmmRJMCMC(x, mInit = 1, mVec = c(1, 2), qnew = 1, burn = 0, niter = 1, verbose = FALSE),
      NA
    ),
    "deprecated"
  )
})

test_that("pgmm_rjmcmc honors zero posterior iterations", {
  set.seed(2026)
  x <- cbind(
    matrix(rnorm(8, mean = -2, sd = 0.2), nrow = 2),
    matrix(rnorm(8, mean = 2, sd = 0.2), nrow = 2)
  )

  fit <- pgmm_rjmcmc(
    x,
    mInit = 2,
    mVec = c(1, 3),
    qnew = 1,
    burn = 0,
    niter = 0,
    verbose = FALSE
  )

  expect_length(fit$ZmatList, 0)
  expect_length(fit$constraintList, 0)
  expect_length(fit$alpha1Vec, 0)
})

test_that("pgmm_rjmcmc validates user-facing inputs", {
  x <- matrix(rnorm(12), nrow = 2)

  expect_error(pgmm_rjmcmc(x, mInit = 0, mVec = c(1, 2), qnew = 1), "mInit")
  expect_error(pgmm_rjmcmc(x, mInit = 1, mVec = c(2, 1), qnew = 1), "mVec")
  expect_error(pgmm_rjmcmc(x, mInit = 1, mVec = c(1, 2), qnew = 0), "qnew")
  expect_error(pgmm_rjmcmc(x, mInit = 1, mVec = c(1, 2), qnew = 1, burn = -1), "burn")
  expect_error(pgmm_rjmcmc(x, mInit = 1, mVec = c(1, 2), qnew = 1, niter = -1), "niter")
  expect_error(pgmm_rjmcmc(x, mInit = 1, mVec = c(1, 2), qnew = 1, constraint = c(0, 0)), "constraint")
  expect_error(pgmm_rjmcmc(x, mInit = 1, mVec = c(1, 2), qnew = 1, Mstep = 2), "Mstep")
  expect_error(pgmm_rjmcmc(x, mInit = 1, mVec = c(1, 2), qnew = 1, Vstep = NA), "Vstep")
  expect_error(pgmm_rjmcmc(x, mInit = 1, mVec = c(1, 2), qnew = 1, SCind = -1), "SCind")
  expect_error(pgmm_rjmcmc(as.data.frame(x), mInit = 1, mVec = c(1, 2), qnew = 1), "X")
})
