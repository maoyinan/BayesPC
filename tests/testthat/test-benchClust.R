test_that("multiplication works", {
  ls_tb <- benchClust(DATASET, "ID", "t", "Record", "Group", out_pc$arr_cluster)
  expect(all(dplyr::between(range(unlist(ls_tb), na.rm = T), 0, 1)), failure_message = "Rand index out of range.")
})
