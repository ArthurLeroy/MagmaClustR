test_that("gradient of gr_clust_multi_GP() works for Squared Exponential kernel", {
  k = seq_len(3)
  m_k <- c("K1" = 0, "K2" = 0, "K3" = 0)

  db <- simu_db(N = 10, common_input = TRUE)
  hp_k <- MagmaClustR:::hp("SE", list_ID = names(m_k))
  hp_i <- MagmaClustR:::hp("SE", list_ID = unique(db$ID))

  old_tau_i_k = ini_tau_i_k(db = db, k = length(k), nstart = 50)
  hp_k[['pi']] = sapply( old_tau_i_k, function(x) x %>% unlist() %>% mean() )

  mu_k_param = MagmaClustR:::e_step_VEM(db, m_k, "SE", "SE", hp_k, hp_i, old_tau_i_k ,0.001)

  db_1 <- db %>% dplyr::filter(.data$ID == 1)
  hp_i <- hp_i %>% dplyr::select(-.data$ID)
  hp_i_1 <- hp_i[1,]

  hp_v = hp_l = hp_i_1
  hp_v$variance = hp_i_1$variance +  10^(-8)
  hp_l$lengthscale = hp_i_1$lengthscale +  10^(-8)

  deriv_v = gr_clust_multi_GP(hp_i_1, db_1, mu_k_param, 'SE', 0)[['variance']]
  deriv_l = gr_clust_multi_GP(hp_i_1, db_1, mu_k_param, 'SE', 0)[['lengthscale']]

  emp_deriv_v = (logL_clust_multi_GP(hp_v, db_1, mu_k_param, 'SE', 0) -
                   logL_clust_multi_GP(hp_i_1, db_1, mu_k_param, 'SE', 0)) / 10^(-8)
  emp_deriv_l = (logL_clust_multi_GP(hp_l, db_1, mu_k_param, 'SE', 0) -
                   logL_clust_multi_GP(hp_i_1, db_1, mu_k_param, 'SE', 0)) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(emp_deriv_v, 4))
  round(deriv_l, 4) %>% expect_equal(round(emp_deriv_l, 4))
})




