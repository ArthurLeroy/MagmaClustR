test_that("gradient of gr_clust_multi_GP() works for Squared Exponential kernel", {
  k = seq_len(3)
  m_k <- c("K1" = 0, "K2" = 0, "K3" = 0)

  db <- simu_db(N = 10, common_input = TRUE)
  hp_k <- MagmaClustR:::hp("SE", list_ID = names(m_k))
  hp_i <- MagmaClustR:::hp("SE", list_ID = unique(db$ID))


  #db <- tibble::tibble('ID' = 1:5,Input= 1:5, Output= 2:6, Covariates = 3:7)
  #hp_k <- tibble::tibble('ID' = names(m_k), 'se_variance' = 1, 'se_lengthscale' = 1)
  #hp_i <- tibble::tibble('ID' = unique(db$ID), 'se_variance' = 1, 'se_lengthscale' = 1)

  old_hp_mixture = ini_hp_mixture(db = db, k = length(k), nstart = 50)
  prop_mixture_1 <- old_hp_mixture %>% dplyr::select(-.data$ID)
  hp_k[['prop_mixture']] = sapply( prop_mixture_1, function(x) x %>% unlist() %>% mean() )

  mu_k_param = MagmaClustR:::e_step_VEM(db, m_k, "SE", "SE", hp_k, hp_i, old_hp_mixture ,0.001)

  db_1 <- db %>% dplyr::filter(.data$ID == 1)
  hp_i_1 <- hp_i %>% dplyr::select(-.data$ID) %>% dplyr::slice(1)

  hp_v = hp_l = hp_i_1
  hp_v$se_variance = hp_i_1$se_variance +  10^(-8)
  hp_l$se_lengthscale = hp_i_1$se_lengthscale +  10^(-8)

  deriv_v = gr_clust_multi_GP(hp_i_1, db_1, mu_k_param, 'SE', 1)[['se_variance']]
  deriv_l = gr_clust_multi_GP(hp_i_1, db_1, mu_k_param, 'SE', 1)[['se_lengthscale']]

  emp_deriv_v = (elbo_clust_multi_GP(hp_v, db_1, mu_k_param, 'SE', 1) -
                   elbo_clust_multi_GP(hp_i_1, db_1, mu_k_param, 'SE', 1)) / 10^(-8)
  emp_deriv_l = (elbo_clust_multi_GP(hp_l, db_1, mu_k_param, 'SE', 1) -
                   elbo_clust_multi_GP(hp_i_1, db_1, mu_k_param, 'SE', 1)) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(as.double(emp_deriv_v), 4))
  round(deriv_l, 4) %>% expect_equal(round(as.double(emp_deriv_l), 4))
})

test_that("gradient of gr_GP_mod_common_hp_k() works", {

  k = seq_len(3)
  m_k <- c("K1" = 0, "K2" = 0, "K3" = 0)

  #db <- simu_db(N = 10, common_input = TRUE)
  #hp_k <- MagmaClustR:::hp("SE", list_ID = names(m_k))
  #hp_i <- MagmaClustR:::hp("SE", list_ID = unique(db$ID))


  db <- tibble::tibble('ID' = 1:5,Input= 1:5, Output= 2:6, Covariates = 3:7)
  hp_k <- tibble::tibble('ID' = names(m_k), 'se_variance' = 1, 'se_lengthscale' = 1)
  hp_i <- tibble::tibble('ID' = unique(db$ID), 'se_variance' = 1, 'se_lengthscale' = 1)

  old_hp_mixture = ini_hp_mixture(db = db, k = length(k), nstart = 50)
  prop_mixture_1 <- old_hp_mixture %>% dplyr::select(-.data$ID)
  hp_k[['prop_mixture']] = sapply( prop_mixture_1, function(x) x %>% unlist() %>% mean() )

  mu_k_param = MagmaClustR:::e_step_VEM(db, m_k, "SE", "SE", hp_k, hp_i, old_hp_mixture ,0.001)

  #db_1 <- db %>% dplyr::filter(.data$ID == 1)
  hp_i_1 <- hp_i %>% dplyr::select(-.data$ID) %>% dplyr::slice(1)

  hp_v = hp_l = hp_i_1
  hp_v$se_variance = hp_i_1$se_variance +  10^(-8)
  hp_l$se_lengthscale = hp_i_1$se_lengthscale +  10^(-8)


  db_1 <- mu_k_param$mean
  new_cov = mu_k_param$cov


  deriv_v = gr_GP_mod_common_hp_k(hp_i_1, db_1, m_k, 'SE', new_cov, 1)[['se_variance']] %>% sum
  deriv_l = gr_GP_mod_common_hp_k(hp_i_1, db_1, m_k, 'SE', new_cov, 1)[['se_lengthscale']] %>% sum

  emp_deriv_v = (elbo_GP_mod_common_hp_k(hp_v, db_1, m_k, 'SE', new_cov, 1) -
                   elbo_GP_mod_common_hp_k(hp_i_1, db_1, m_k, 'SE', new_cov, 1)) / 10^(-8)

  emp_deriv_l = (elbo_GP_mod_common_hp_k(hp_l, db_1, m_k, 'SE', new_cov, 1) -
                   elbo_GP_mod_common_hp_k(hp_i_1, db_1, m_k, 'SE', new_cov, 1)) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(emp_deriv_v, 4))
  round(deriv_l, 4) %>% expect_equal(round(emp_deriv_l, 4))
})

test_that("gradient of gr_clust_multi_GP_common_hp_i() works", {
  k = seq_len(3)
  m_k <- c("K1" = 0, "K2" = 0, "K3" = 0)

  #db <- simu_db(N = 10, common_input = TRUE)
  #hp_k <- MagmaClustR:::hp("SE", list_ID = names(m_k))
  #hp_i <- MagmaClustR:::hp("SE", list_ID = unique(db$ID))


  db <- tibble::tibble('ID' = 1:5,Input= 1:5, Output= 2:6, Covariates = 3:7)
  hp_k <- tibble::tibble('ID' = names(m_k), 'se_variance' = 1, 'se_lengthscale' = 1)
  hp_i <- tibble::tibble('ID' = unique(db$ID), 'se_variance' = 1, 'se_lengthscale' = 1)

  old_hp_mixture = ini_hp_mixture(db = db, k = length(k), nstart = 50)
  prop_mixture_1 <- old_hp_mixture %>% dplyr::select(-.data$ID)
  hp_k[['prop_mixture']] = sapply( prop_mixture_1, function(x) x %>% unlist() %>% mean() )

  mu_k_param = MagmaClustR:::e_step_VEM(db, m_k, "SE", "SE", hp_k, hp_i, old_hp_mixture ,0.001)

  db_1 <- db %>% dplyr::filter(.data$ID == 1)
  hp_i_1 <- hp_i %>% dplyr::select(-.data$ID) %>% dplyr::slice(1)

  hp_v = hp_l = hp_i_1
  hp_v$se_variance = hp_i_1$se_variance +  10^(-8)
  hp_l$se_lengthscale = hp_i_1$se_lengthscale +  10^(-8)



  deriv_v = gr_clust_multi_GP_common_hp_i(hp_i_1, db_1, mu_k_param, 'SE', 1)[['se_variance']]
  deriv_l = gr_clust_multi_GP_common_hp_i(hp_i_1, db_1, mu_k_param, 'SE', 1)[['se_lengthscale']]

  emp_deriv_v = (elbo_clust_multi_GP_common_hp_i(hp_v, db_1, mu_k_param, 'SE', 1) -
                   elbo_clust_multi_GP_common_hp_i(hp_i_1, db_1, mu_k_param, 'SE', 1)) / 10^(-8)
  emp_deriv_l = (elbo_clust_multi_GP_common_hp_i(hp_l, db_1, mu_k_param, 'SE', 1) -
                   elbo_clust_multi_GP_common_hp_i(hp_i_1, db_1, mu_k_param, 'SE', 1)) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(as.double(emp_deriv_v), 4))
  round(deriv_l, 4) %>% expect_equal(round(as.double(emp_deriv_l), 4))
})


