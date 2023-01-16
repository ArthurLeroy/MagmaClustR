test_that("simu_db() returns a dataframe", {
  simu_db() %>% expect_s3_class('data.frame')
})

