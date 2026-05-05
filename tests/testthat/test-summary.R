test_that("summarizePgmmRJMCMC summarizes posterior allocation and models", {
  pgmm_res <- list(
    ZmatList = list(c(1, 1, 2), c(1, 2, 2), c(1, 1, 2)),
    constraintList = list(c(0, 0, 0), c(1, 0, 0), c(1, 0, 0))
  )

  summary <- summarizePgmmRJMCMC(pgmm_res, trueCluster = c(1, 1, 2))

  expect_equal(summary$Zalloc, c(1, 1, 2))
  expect_equal(as.integer(summary$nCluster["2"]), 3L)
  expect_equal(as.integer(summary$nConstraint["CUU"]), 2L)
  expect_equal(as.integer(summary$nConstraint["UUU"]), 1L)
  expect_equal(summary$ari, 1)
})

test_that("misspelled summary function remains a compatibility alias", {
  pgmm_res <- list(
    ZmatList = list(c(1, 2), c(1, 2)),
    constraintList = list(c(0, 0, 0), c(0, 0, 0))
  )

  expect_equal(
    summerizePgmmRJMCMC(pgmm_res),
    summarizePgmmRJMCMC(pgmm_res)
  )
})
