test_that("simu_db() returns a dataframe", {
  simu_db() %>% inherits('data.frame')
})

