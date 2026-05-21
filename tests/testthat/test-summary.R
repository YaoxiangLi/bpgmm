test_that("summarize_pgmm_rjmcmc summarizes posterior allocation and models", {
  pgmm_res <- list(
    ZmatList = list(c(1, 1, 2), c(1, 2, 2), c(1, 1, 2)),
    constraintList = list(c(0, 0, 0), c(1, 0, 0), c(1, 0, 0))
  )

  summary <- summarize_pgmm_rjmcmc(pgmm_res, trueCluster = c(1, 1, 2))

  expect_equal(summary$Zalloc, c(1, 1, 2))
  expect_equal(as.integer(summary$nCluster["2"]), 3L)
  expect_equal(as.integer(summary$nConstraint["CUU"]), 2L)
  expect_equal(as.integer(summary$nConstraint["UUU"]), 1L)
  expect_equal(summary$ari, 1)
})

test_that("old summary function names are deprecated compatibility aliases", {
  pgmm_res <- list(
    ZmatList = list(c(1, 2), c(1, 2)),
    constraintList = list(c(0, 0, 0), c(0, 0, 0))
  )

  expect_warning(
    old_summary <- summarizePgmmRJMCMC(pgmm_res),
    "deprecated"
  )
  expect_warning(
    misspelled_summary <- summerizePgmmRJMCMC(pgmm_res),
    "deprecated"
  )
  expect_equal(old_summary, summarize_pgmm_rjmcmc(pgmm_res))
  expect_equal(misspelled_summary, summarize_pgmm_rjmcmc(pgmm_res))
})

test_that("internal allocation summarizer keeps compatibility alias", {
  z_samples <- list(c(1, 1, 2), c(1, 2, 2), c(1, 1, 2))

  expect_equal(bpgmm:::summarizeZ(z_samples), c(1, 1, 2))
  expect_equal(bpgmm:::sumerizeZ(z_samples), bpgmm:::summarizeZ(z_samples))
})

test_that("summarize_pgmm_rjmcmc validates result structure", {
  expect_error(summarize_pgmm_rjmcmc(list()), "ZmatList")
  expect_error(
    summarize_pgmm_rjmcmc(list(ZmatList = list(), constraintList = list())),
    "ZmatList"
  )
  expect_error(
    summarize_pgmm_rjmcmc(list(ZmatList = list(c(1, 2)), constraintList = list())),
    "constraintList"
  )
})
