test_that("Parameter calculation works", {
  ls_par <- postMean(df_of_draws, paste0("Z", 1:10), "ID", DATASET)
  nSub <- length(unique(DATASET[, "ID"]))

  expect_equal(nrow(ls_par$bMatrix), nSub)
  expect_equal(ncol(ls_par$bMatrix), nrow(df_of_draws))
  expect_equal(nrow(ls_par$G), 10)
})
