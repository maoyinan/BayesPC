test_that("LMM using stan works", {
  expect_equal(
    nrow(df_of_draws), 10
  )
})
