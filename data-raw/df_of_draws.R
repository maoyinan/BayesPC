## code to prepare `df_of_draws` dataset goes here

data(DATASET)
df_of_draws <- modelStan("Record", paste0("Z", 1:10), "ID", DATASET)
df_of_draws <- df_of_draws[1:10, ]
usethis::use_data(df_of_draws, overwrite = TRUE)
