test_that("Prediction calculation works", {
  ls_par <- postMean(df_of_draws,paste0('Z',1:10), 'ID', DATASET)
  ls_idxA <- list(
   seq(10),
   1:4,
   5:7,
   8:10
  )
  mat_fitted <- pcFit('ID', ls_par, DATASET, ls_idxA)

  expect_equal(length(ls_idxA), ncol(mat_fitted))
})
