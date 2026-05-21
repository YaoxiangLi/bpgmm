test_that("Rcpp get_z_mat_cpp matches the R implementation", {
  labels <- c(1, 2, 1, 3)

  expect_equal(
    bpgmm:::get_z_mat_cpp(labels, m = 3, n = length(labels)),
    bpgmm:::get_z_mat_r(labels, m = 3, n = length(labels))
  )
})

test_that("Rcpp get_z_mat_cpp rejects labels outside 1:m", {
  expect_error(
    bpgmm:::get_z_mat_cpp(c(1, 0, 2), m = 2, n = 3),
    "cluster labels must be integers in 1:m"
  )
  expect_error(
    bpgmm:::get_z_mat_cpp(c(1, 3, 2), m = 2, n = 3),
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

test_that("Rcpp dmvnrm_arma validates dimensions and covariance", {
  x <- matrix(c(1, 2), nrow = 1)

  expect_error(
    bpgmm:::dmvnrm_arma(x, mean = 0, sigma = diag(2), logd = TRUE),
    "mean length"
  )
  expect_error(
    bpgmm:::dmvnrm_arma(x, mean = c(0, 0), sigma = diag(3), logd = TRUE),
    "sigma"
  )
  expect_error(
    bpgmm:::dmvnrm_arma(
      x,
      mean = c(0, 0),
      sigma = matrix(c(1, 2, 2, 1), nrow = 2),
      logd = TRUE
    ),
    "positive definite"
  )
})

test_that("Rcpp calculate_ratio_cpp is stable for large log probabilities", {
  expect_equal(
    bpgmm:::calculate_ratio_cpp(1000, c(1000, 1001, 1002)),
    exp(1000 - 1002) / sum(exp(c(1000, 1001, 1002) - 1002)),
    tolerance = 1e-12
  )
})

test_that("Rcpp calculate_ratio_cpp validates finite log inputs", {
  expect_error(
    bpgmm:::calculate_ratio_cpp(0, numeric()),
    "logNume"
  )
  expect_error(
    bpgmm:::calculate_ratio_cpp(Inf, c(0, 1)),
    "finite"
  )
  expect_error(
    bpgmm:::calculate_ratio_cpp(0, c(0, NA)),
    "finite"
  )
})

test_that("Rcpp update_post_z_cpp uses log mixture weights", {
  n <- 2000
  theta <- new(
    "ThetaYList",
    tao = c(0.99, 0.01),
    psy = list(matrix(1), matrix(1)),
    M = list(0, 0),
    lambda = list(matrix(0), matrix(0)),
    Y = list(matrix(0, nrow = 1, ncol = n), matrix(0, nrow = 1, ncol = n))
  )

  set.seed(22)
  labels <- bpgmm:::update_post_z_cpp(matrix(0, nrow = 1, ncol = n), m = 2, n = n, theta)

  expect_gt(mean(labels == 1), 0.95)
})

test_that("Rcpp update_post_z_cpp validates dimensions and mixture weights", {
  theta <- new(
    "ThetaYList",
    tao = c(1, 0),
    psy = list(matrix(1), matrix(1)),
    M = list(0, 0),
    lambda = list(matrix(0), matrix(0)),
    Y = list(matrix(0, nrow = 1, ncol = 1), matrix(0, nrow = 1, ncol = 1))
  )

  expect_error(
    bpgmm:::update_post_z_cpp(matrix(0, nrow = 1, ncol = 1), m = 2, n = 1, theta),
    "tao"
  )
  expect_error(
    bpgmm:::update_post_z_cpp(matrix(0, nrow = 1, ncol = 2), m = 2, n = 1, theta),
    "n"
  )
})
