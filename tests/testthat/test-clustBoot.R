test_that("multiplication works", {
  ls_par <- postMean(df_of_draws, paste0("Z", 1:10), "ID", DATASET)
  ls_idxA <- list(
    seq(10),
    1:4,
    5:7,
    8:10
  )
  mat_fitted <- pcFit("ID", ls_par, DATASET, ls_idxA)
  out_boot1 <- clustBoot(paste0("Z", 1:10), "ID","t", mat_fitted,DATASET, 10, 10, 1)
  out_boot2 <- clustBoot(paste0("Z", 1:10), "ID","t", mat_fitted,DATASET, 10, 10, 1)
  expect_equal(out_boot1, out_boot2)
})
