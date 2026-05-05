test_that("documented public API remains exported", {
  exports <- getNamespaceExports("bpgmm")

  expect_true(all(c(
    "pgmmRJMCMC",
    "constraint_to_model",
    "model_to_constraint",
    "summarizePgmmRJMCMC",
    "summerizePgmmRJMCMC"
  ) %in% exports))
})
