test_that("gradient of gr_clust_multi_GP() works for Squared Exponential kernel", {
  k = seq_len(3)
  m_k <- c("K1" = 0, "K2" = 0, "K3" = 0)

  db <- simu_db(N = 10, common_input = TRUE)
  #db <- tibble::tibble('ID' = 1:5,Input= 1:5, Output= 2:6, Covariates = 3:7)
  #hp_k <- MagmaClustR:::hp("SE", list_ID = names(m_k))
  hp_k <- tibble::tibble('ID' = names(m_k), 'variance' = 1, 'lengthscale' = 1)
  hp_i <- tibble::tibble('ID' = unique(db$ID), 'variance' = 1, 'lengthscale' = 1)
  #hp_i <- MagmaClustR:::hp("SE", list_ID = unique(db$ID))

  old_tau_i_k = ini_tau_i_k(db = db, k = length(k), nstart = 50)
  hp_k[['pi']] = sapply( old_tau_i_k, function(x) x %>% unlist() %>% mean() )

  mu_k_param = MagmaClustR:::e_step_VEM(db, m_k, "SE", "SE", hp_k, hp_i, old_tau_i_k ,0.001)

  db_1 <- db %>% dplyr::filter(.data$ID == 1)
  hp_i_1 <- hp_i %>% dplyr::select(-.data$ID) %>% dplyr::slice(1)
  #hp_i_1 <- hp_i[1,]

  hp_v = hp_l = hp_i_1
  hp_v$variance = hp_i_1$variance +  10^(-8)
  hp_l$lengthscale = hp_i_1$lengthscale +  10^(-8)

  deriv_v = gr_clust_multi_GP(hp_i_1, db_1, mu_k_param, 'SE', 1)[['variance']]
  deriv_l = gr_clust_multi_GP(hp_i_1, db_1, mu_k_param, 'SE', 1)[['lengthscale']]

  emp_deriv_v = (logL_clust_multi_GP(hp_v, db_1, mu_k_param, 'SE', 1) -
                   logL_clust_multi_GP(hp_i_1, db_1, mu_k_param, 'SE', 1)) / 10^(-8)
  emp_deriv_l = (logL_clust_multi_GP(hp_l, db_1, mu_k_param, 'SE', 1) -
                   logL_clust_multi_GP(hp_i_1, db_1, mu_k_param, 'SE', 1)) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(emp_deriv_v, 4))
  round(deriv_l, 4) %>% expect_equal(round(emp_deriv_l, 4))
})

test_that("gradient of gr_GP_mod_common_hp_k() works", {
  # db = tibble::tibble(ID = rep(1:5, each = 4),
  #                     Output = 1:20,
  #                     Input = 2:21,
  #                     Covariate = c(1:10, 23, 77, 1:8))
  # mean <- tibble::tibble(Input = unique(db$Input), Output = 0)
  # hp = tibble::tibble(variance = 1, lengthscale = 0.5)

  k = seq_len(3)
  m_k <- c("K1" = 0, "K2" = 0, "K3" = 0)

  db <- simu_db(N = 10, common_input = TRUE)
  hp_k <- tibble::tibble('ID' = names(m_k), 'variance' = 1, 'lengthscale' = 0.5)
  hp_i <- tibble::tibble('ID' = unique(db$ID), 'variance' = 1, 'lengthscale' = 0.5)

  old_tau_i_k = ini_tau_i_k(db = db, k = length(k), nstart = 50)
  hp_k[['pi']] = sapply( old_tau_i_k, function(x) x %>% unlist() %>% mean() )

  mu_k_param = MagmaClustR:::e_step_VEM(db, m_k, "SE", "SE", hp_k, hp_i, old_tau_i_k ,0.001)
  hp_i_1 <- hp_i %>% dplyr::slice(1)


  db_1 <- mu_k_param$mean
  new_cov = mu_k_param$cov

  hp_v = tibble::tibble(variance = 1 + 10^(-8), lengthscale = 0.5)
  hp_l = tibble::tibble(variance = 1, lengthscale = 0.5 + 10^(-8))

  deriv_v = gr_GP_mod_common_hp_k(hp_i_1, db_1, m_k, 'SE', new_cov, 1)['variance',] %>% sum
  deriv_l = gr_GP_mod_common_hp_k(hp_i_1, db_1, m_k, 'SE', new_cov, 1)['lengthscale',] %>% sum

  emp_deriv_v = (logL_GP_mod_common_hp_k(hp_v, db_1, m_k, 'SE', new_cov, 1) -
                   logL_GP_mod_common_hp_k(hp_i_1, db_1, m_k, 'SE', new_cov, 1)) / 10^(-8)

  emp_deriv_l = (logL_GP_mod_common_hp_k(hp_l, db_1, m_k, 'SE', new_cov, 1) -
                   logL_GP_mod_common_hp_k(hp_i_1, db_1, m_k, 'SE', new_cov, 1)) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(emp_deriv_v, 4))
  round(deriv_l, 4) %>% expect_equal(round(emp_deriv_l, 4))
})

test_that("gradient of gr_clust_multi_GP_common_hp_i() works", {
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

  deriv_v = gr_clust_multi_GP_common_hp_i(hp_i_1, db_1, mu_k_param, 'SE', 1)[['variance']]
  deriv_l = gr_clust_multi_GP_common_hp_i(hp_i_1, db_1, mu_k_param, 'SE', 1)[['lengthscale']]

  emp_deriv_v = (logL_clust_multi_GP_common_hp_i(hp_v, db_1, mu_k_param, 'SE', 1) -
                   logL_clust_multi_GP_common_hp_i(hp_i_1, db_1, mu_k_param, 'SE', 1)) / 10^(-8)
  emp_deriv_l = (logL_clust_multi_GP_common_hp_i(hp_l, db_1, mu_k_param, 'SE', 1) -
                   logL_clust_multi_GP_common_hp_i(hp_i_1, db_1, mu_k_param, 'SE', 1)) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(emp_deriv_v, 4))
  round(deriv_l, 4) %>% expect_equal(round(emp_deriv_l, 4))
})


