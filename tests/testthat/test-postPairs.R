test_that("Probability matrix valued between 0 and 1", {
  ls_idxA <- list(
    seq(10),
    1:4,
    5:7,
    8:10
  )
  ls_prob <- postPairs(df_of_draws,
    x_var = paste0("Z", 1:10), id_var = "ID", dat = DATASET,
    ls_idxA, nIter = 10, nDraw = 2, nClusters = 4, regQ = 1e-6, seed = 1
  )$ls_prob

  expect(all(dplyr::between(range(unlist(ls_prob), na.rm = T), 0, 1)), failure_message = "Probability value out of range.")
})
