test_that("random seed works", {
  ls_par <- postMean(df_of_draws, paste0("Z", 1:10), "ID", DATASET)
  ls_idxA <- list(
    seq(10),
    1:4,
    5:7,
    8:10
  )
  mat_fitted <- pcFit("ID", ls_par, DATASET, ls_idxA)
  out_KL1 <- clustKL("ID", ls_par, DATASET, ls_idxA, 10, .1, 1e-6, 1)
  out_KL2 <- clustKL("ID", ls_par, DATASET, ls_idxA, 10, .1, 1e-6, 1)
  out_KL3 <- clustKL("ID", ls_par, DATASET, ls_idxA, 10, .1, 1e-6, 3)

  expect_equal(out_KL1, out_KL2)
})
