test_that("list_kern_to_cov() and list_kern_to_inv() are inverse of the other", {

  db = simu_db(M = 2, N = 3)
  hp = hp('SE', unique(db$ID))

  list_cov = list_kern_to_cov(db, 'SE', hp)
  list_inv = list_kern_to_inv(db, 'SE', hp)
  for(i in unique(db$ID)){
    list_cov[[i]] %>% expect_equal(list_inv[[i]] %>% solve())
  }

})
