test_that("Rcpp get_Z_mat matches the R implementation", {
  labels <- c(1, 2, 1, 3)

  expect_equal(
    bpgmm:::get_Z_mat(labels, m = 3, n = length(labels)),
    bpgmm:::getZmat(labels, m = 3, n = length(labels))
  )
})

test_that("Rcpp get_Z_mat rejects labels outside 1:m", {
  expect_error(
    bpgmm:::get_Z_mat(c(1, 0, 2), m = 2, n = 3),
    "cluster labels must be integers in 1:m"
  )
  expect_error(
    bpgmm:::get_Z_mat(c(1, 3, 2), m = 2, n = 3),
    "cluster labels must be integers in 1:m"
  )
})

test_that("Rcpp dmvnrm_arma matches mvtnorm::dmvnorm", {
  x <- rbind(c(0, 0), c(1, -1))
  mean <- c(0.25, -0.5)
  sigma <- matrix(c(2, 0.4, 0.4, 1.5), nrow = 2)

  expect_equal(
    as.numeric(bpgmm:::dmvnrm_arma(x, mean, sigma, TRUE)),
    as.numeric(mvtnorm::dmvnorm(x, mean = mean, sigma = sigma, log = TRUE)),
    tolerance = 1e-10
  )
})

test_that("Rcpp calculate_Ratio is stable for large log probabilities", {
  expect_equal(
    bpgmm:::calculate_Ratio(1000, c(1000, 1001, 1002)),
    exp(1000 - 1002) / sum(exp(c(1000, 1001, 1002) - 1002)),
    tolerance = 1e-12
  )
})
