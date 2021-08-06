## code to prepare `ls_prob` dataset goes here

data(DATASET)
df_of_draws <- modelStan("Record", paste0("Z", 1:10), "ID", DATASET, seed = 1)
ls_idxA <- list(
  seq(10),
  1:4,
  5:7,
  8:10
)
out_pc <- postPairs(df_of_draws,
  x_var = paste0("Z", 1:10), id_var = "ID", dat = DATASET,
  ls_idxA, nIter = 10, nDraw = 1000, nClusters = 4, regQ = 1e-6, seed = 1
)
usethis::use_data(out_pc, overwrite = TRUE)
