test_that("calculateRatio is stable for log-scale probabilities", {
  expect_equal(
    bpgmm:::calculateRatio(log(2), log(c(2, 3, 5))),
    0.2,
    tolerance = 1e-12
  )
})

test_that("getZmat converts one-dimensional labels to indicator rows", {
  expect_equal(
    bpgmm:::getZmat(c(1, 2, 1, 3), m = 3, n = 4),
    rbind(c(1, 0, 1, 0), c(0, 1, 0, 0), c(0, 0, 0, 1))
  )
})
