########################## gr_clust_multi_GP() #####################

test_that("gradient of gr_clust_multi_GP() works for
          the Squared Exponential kernel", {
  k <- seq_len(3)
  m_k <- list(
    "K1" = rep(0, 5),
    "K2" = rep(0, 5),
    "K3" = rep(0, 5)
  )

  # db <- simu_db(N = 10, common_input = TRUE)
  # hp_k <- MagmaClustR:::hp("SE", list_ID = names(m_k))
  # hp_i <- MagmaClustR:::hp("SE", list_ID = unique(db$ID))


  db <- tibble::tibble("ID" = 1:5, Input = 1:5, Output = 2:6, Covariates = 3:7)
  hp_k <- tibble::tibble(
    "ID" = names(m_k),
    "se_variance" = 1,
    "se_lengthscale" = 1
  )
  hp_i <- tibble::tibble(
    "ID" = unique(db$ID),
    "se_variance" = 1,
    "se_lengthscale" = 1
  )

  old_hp_mixture <- ini_mixture(data = db, k = length(k), nstart = 50)
  prop_mixture_1 <- old_hp_mixture %>% dplyr::select(-.data$ID)
  hp_k[["prop_mixture"]] <- sapply(prop_mixture_1, function(x) {
    x %>%
      unlist() %>%
      mean()
  })

  mu_k_param <- MagmaClustR:::ve_step(
    db, m_k, "SE", "SE",
    hp_k, hp_i, old_hp_mixture, 2, 0.001
  )

  db_1 <- db %>% dplyr::filter(.data$ID == 1)
  hp_i_1 <- hp_i %>%
    dplyr::select(-.data$ID) %>%
    dplyr::slice(1)

  hp_v <- hp_l <- hp_i_1
  hp_v$se_variance <- hp_i_1$se_variance + 10^(-8)
  hp_l$se_lengthscale <- hp_i_1$se_lengthscale + 10^(-8)

  deriv_v <- gr_clust_multi_GP(
    hp_i_1, db_1,
    mu_k_param, "SE", 1
  )[["se_variance"]]
  deriv_l <- gr_clust_multi_GP(
    hp_i_1, db_1,
    mu_k_param, "SE", 1
  )[["se_lengthscale"]]

  emp_deriv_v <- (elbo_clust_multi_GP(hp_v, db_1, mu_k_param, "SE", 1) -
    elbo_clust_multi_GP(hp_i_1, db_1, mu_k_param, "SE", 1)) / 10^(-8)
  emp_deriv_l <- (elbo_clust_multi_GP(hp_l, db_1, mu_k_param, "SE", 1) -
    elbo_clust_multi_GP(hp_i_1, db_1, mu_k_param, "SE", 1)) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(as.double(emp_deriv_v), 4))
  round(deriv_l, 4) %>% expect_equal(round(as.double(emp_deriv_l), 3))
})

test_that("gradient of gr_clust_multi_GP() works for the Linear kernel", {
  k <- seq_len(3)
  m_k <- list("K1" = rep(0, 5), "K2" = rep(0, 5), "K3" = rep(0, 5))

  # db <- simu_db(N = 10, common_input = TRUE)
  # hp_k <- MagmaClustR:::hp("LIN", list_ID = names(m_k))
  # hp_i <- MagmaClustR:::hp("LIN", list_ID = unique(db$ID))


  db <- tibble::tibble("ID" = 1:5, Input = 1:5, Output = 2:6, Covariates = 3:7)
  hp_k <- tibble::tibble("ID" = names(m_k), "lin_slope" = 1, "lin_offset" = 1)
  hp_i <- tibble::tibble("ID" = unique(db$ID), "lin_slope" = 1, "lin_offset" = 1)

  old_hp_mixture <- ini_mixture(data = db, k = length(k), nstart = 50)
  prop_mixture_1 <- old_hp_mixture %>% dplyr::select(-.data$ID)
  hp_k[["prop_mixture"]] <- sapply(
    prop_mixture_1,
    function(x) {
      x %>%
        unlist() %>%
        mean()
    }
  )

  mu_k_param <- MagmaClustR:::ve_step(
    db, m_k, "LIN", "LIN",
    hp_k, hp_i, old_hp_mixture, 2, 0.001
  )

  db_1 <- db %>% dplyr::filter(.data$ID == 1)
  hp_i_1 <- hp_i %>%
    dplyr::select(-.data$ID) %>%
    dplyr::slice(1)

  hp_v <- hp_l <- hp_i_1
  hp_v$lin_slope <- hp_i_1$lin_slope + 10^(-10)
  hp_l$lin_offset <- hp_i_1$lin_offset + 10^(-10)

  deriv_v <- gr_clust_multi_GP(hp_i_1, db_1, mu_k_param, "LIN", 1)[["lin_slope"]]
  deriv_l <- gr_clust_multi_GP(hp_i_1, db_1, mu_k_param, "LIN", 1)[["lin_offset"]]

  emp_deriv_v <- (elbo_clust_multi_GP(hp_v, db_1, mu_k_param, "LIN", 1) -
    elbo_clust_multi_GP(hp_i_1, db_1, mu_k_param, "LIN", 1)) / 10^(-10)
  emp_deriv_l <- (elbo_clust_multi_GP(hp_l, db_1, mu_k_param, "LIN", 1) -
    elbo_clust_multi_GP(hp_i_1, db_1, mu_k_param, "LIN", 1)) / 10^(-10)

  round(deriv_v, 3) %>% expect_equal(round(as.double(emp_deriv_v), 3))
  round(deriv_l, 3) %>% expect_equal(round(as.double(emp_deriv_l), 3))
})

test_that("gradient of gr_clust_multi_GP() works for the Periodic kernel", {
  k <- seq_len(3)
  m_k <- list("K1" = rep(0, 5), "K2" = rep(0, 5), "K3" = rep(0, 5))

  # db <- simu_db(N = 10, common_input = TRUE)
  # hp_k <- MagmaClustR:::hp("PERIO", list_ID = names(m_k))
  # hp_i <- MagmaClustR:::hp("PERIO", list_ID = unique(db$ID))


  db <- tibble::tibble("ID" = 1:5, Input = 1:5, Output = 2:6, Covariates = 3:7)
  hp_k <- tibble::tibble(
    "ID" = names(m_k),
    "perio_variance" = 1,
    "perio_lengthscale" = 1,
    "period" = 1
  )
  hp_i <- tibble::tibble(
    "ID" = unique(db$ID),
    "perio_variance" = 1,
    "perio_lengthscale" = 1,
    "period" = 1
  )

  old_hp_mixture <- ini_mixture(data = db, k = length(k), nstart = 50)
  prop_mixture_1 <- old_hp_mixture %>% dplyr::select(-.data$ID)
  hp_k[["prop_mixture"]] <- sapply(prop_mixture_1, function(x) {
    x %>%
      unlist() %>%
      mean()
  })

  mu_k_param <- MagmaClustR:::ve_step(
    db, m_k, "PERIO", "PERIO",
    hp_k, hp_i, old_hp_mixture, 2, 0.001
  )

  db_1 <- db %>% dplyr::filter(.data$ID == 1)
  hp_i_1 <- hp_i %>%
    dplyr::select(-.data$ID) %>%
    dplyr::slice(1)

  hp_v <- hp_l <- hp_p <- hp_i_1
  hp_v$perio_variance <- hp_i_1$perio_variance + 10^(-10)
  hp_l$perio_lengthscale <- hp_i_1$perio_lengthscale + 10^(-10)
  hp_p$period <- hp_i_1$period + 10^(-10)


  deriv_v <- gr_clust_multi_GP(
    hp_i_1, db_1,
    mu_k_param, "PERIO", 1
  )[["perio_variance"]]
  deriv_l <- gr_clust_multi_GP(
    hp_i_1, db_1,
    mu_k_param, "PERIO", 1
  )[["perio_lengthscale"]]
  deriv_p <- gr_clust_multi_GP(
    hp_i_1, db_1,
    mu_k_param, "PERIO", 1
  )[["period"]]


  emp_deriv_v <- (elbo_clust_multi_GP(hp_v, db_1, mu_k_param, "PERIO", 1) -
    elbo_clust_multi_GP(hp_i_1, db_1, mu_k_param, "PERIO", 1)) / 10^(-10)
  emp_deriv_l <- (elbo_clust_multi_GP(hp_l, db_1, mu_k_param, "PERIO", 1) -
    elbo_clust_multi_GP(hp_i_1, db_1, mu_k_param, "PERIO", 1)) / 10^(-10)
  emp_deriv_p <- (elbo_clust_multi_GP(hp_p, db_1, mu_k_param, "PERIO", 1) -
    elbo_clust_multi_GP(hp_i_1, db_1, mu_k_param, "PERIO", 1)) / 10^(-10)

  round(deriv_v, 4) %>% expect_equal(round(as.double(emp_deriv_v), 4))
  round(deriv_l, 4) %>% expect_equal(round(as.double(emp_deriv_l), 4))
  round(deriv_p, 4) %>% expect_equal(round(as.double(emp_deriv_p), 4))
})

test_that("gradient of gr_clust_multi_GP() works
          for the Rational Quadratic kernel", {
  k <- seq_len(3)
  m_k <- list("K1" = rep(0, 5), "K2" = rep(0, 5), "K3" = rep(0, 5))

  # db <- simu_db(N = 10, common_input = TRUE)
  # hp_k <- MagmaClustR:::hp("RQ", list_ID = names(m_k))
  # hp_i <- MagmaClustR:::hp("RQ", list_ID = unique(db$ID))


  db <- tibble::tibble("ID" = 1:5, Input = 1:5, Output = 2:6, Covariates = 3:7)
  hp_k <- tibble::tibble(
    "ID" = names(m_k),
    "rq_variance" = 1,
    "rq_lengthscale" = 1,
    "rq_scale" = 1
  )
  hp_i <- tibble::tibble(
    "ID" = unique(db$ID),
    "rq_variance" = 1,
    "rq_lengthscale" = 1,
    "rq_scale" = 1
  )

  old_hp_mixture <- ini_mixture(data = db, k = length(k), nstart = 50)
  prop_mixture_1 <- old_hp_mixture %>% dplyr::select(-.data$ID)
  hp_k[["prop_mixture"]] <- sapply(
    prop_mixture_1,
    function(x) {
      x %>%
        unlist() %>%
        mean()
    }
  )

  mu_k_param <- MagmaClustR:::ve_step(
    db, m_k, "RQ", "RQ",
    hp_k, hp_i, old_hp_mixture, 2, 0.001
  )

  db_1 <- db %>% dplyr::filter(.data$ID == 1)
  hp_i_1 <- hp_i %>%
    dplyr::select(-.data$ID) %>%
    dplyr::slice(1)

  hp_v <- hp_l <- hp_p <- hp_i_1
  hp_v$rq_variance <- hp_i_1$rq_variance + 10^(-8)
  hp_l$rq_lengthscale <- hp_i_1$rq_lengthscale + 10^(-8)
  hp_p$rq_scale <- hp_i_1$rq_scale + 10^(-8)


  deriv_v <- gr_clust_multi_GP(
    hp_i_1, db_1,
    mu_k_param, "RQ", 1
  )[["rq_variance"]]
  deriv_l <- gr_clust_multi_GP(
    hp_i_1, db_1,
    mu_k_param, "RQ", 1
  )[["rq_lengthscale"]]
  deriv_p <- gr_clust_multi_GP(
    hp_i_1, db_1,
    mu_k_param, "RQ", 1
  )[["rq_scale"]]


  emp_deriv_v <- (elbo_clust_multi_GP(hp_v, db_1, mu_k_param, "RQ", 1) -
    elbo_clust_multi_GP(hp_i_1, db_1, mu_k_param, "RQ", 1)) / 10^(-8)
  emp_deriv_l <- (elbo_clust_multi_GP(hp_l, db_1, mu_k_param, "RQ", 1) -
    elbo_clust_multi_GP(hp_i_1, db_1, mu_k_param, "RQ", 1)) / 10^(-8)
  emp_deriv_p <- (elbo_clust_multi_GP(hp_p, db_1, mu_k_param, "RQ", 1) -
    elbo_clust_multi_GP(hp_i_1, db_1, mu_k_param, "RQ", 1)) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(as.double(emp_deriv_v), 4))
  round(deriv_l, 4) %>% expect_equal(round(as.double(emp_deriv_l), 4))
  round(deriv_p, 4) %>% expect_equal(round(as.double(emp_deriv_p), 4))
})

########################## gr_GP_mod_common_hp_k() #####################

test_that("gradient of gr_GP_mod_common_hp_k()
          works for the Squared Exponential kernel", {
  k <- seq_len(3)
  m_k <- list("K1" = rep(0, 5), "K2" = rep(0, 5), "K3" = rep(0, 5))

  # db <- simu_db(N = 10, common_input = TRUE)
  # hp_k <- MagmaClustR:::hp("SE", list_ID = names(m_k))
  # hp_i <- MagmaClustR:::hp("SE", list_ID = unique(db$ID))


  db <- tibble::tibble("ID" = 1:5, Input = 1:5, Output = 2:6, Covariates = 3:7)
  hp_k <- tibble::tibble(
    "ID" = names(m_k),
    "se_variance" = 1,
    "se_lengthscale" = 1
  )
  hp_i <- tibble::tibble(
    "ID" = unique(db$ID),
    "se_variance" = 1,
    "se_lengthscale" = 1
  )

  old_hp_mixture <- ini_mixture(data = db, k = length(k), nstart = 50)
  prop_mixture_1 <- old_hp_mixture %>% dplyr::select(-.data$ID)
  hp_k[["prop_mixture"]] <- sapply(
    prop_mixture_1,
    function(x) {
      x %>%
        unlist() %>%
        mean()
    }
  )

  mu_k_param <- MagmaClustR:::ve_step(
    db, m_k, "SE", "SE",
    hp_k, hp_i, old_hp_mixture, 2, 0.001
  )

  hp_i_1 <- hp_i %>%
    dplyr::select(-.data$ID) %>%
    dplyr::slice(1)

  hp_v <- hp_l <- hp_i_1
  hp_v$se_variance <- hp_i_1$se_variance + 10^(-8)
  hp_l$se_lengthscale <- hp_i_1$se_lengthscale + 10^(-8)


  db_1 <- mu_k_param$mean
  new_cov <- mu_k_param$cov


  deriv_v <- gr_GP_mod_common_hp_k(
    hp_i_1, db_1, m_k,
    "SE", new_cov, 1
  )[["se_variance"]] %>% sum()
  deriv_l <- gr_GP_mod_common_hp_k(
    hp_i_1, db_1, m_k,
    "SE", new_cov, 1
  )[["se_lengthscale"]] %>% sum()

  emp_deriv_v <- (elbo_GP_mod_common_hp_k(hp_v, db_1, m_k, "SE", new_cov, 1) -
    elbo_GP_mod_common_hp_k(hp_i_1, db_1, m_k, "SE", new_cov, 1)) / 10^(-8)

  emp_deriv_l <- (elbo_GP_mod_common_hp_k(hp_l, db_1, m_k, "SE", new_cov, 1) -
    elbo_GP_mod_common_hp_k(hp_i_1, db_1, m_k, "SE", new_cov, 1)) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(emp_deriv_v, 4))
  round(deriv_l, 4) %>% expect_equal(round(emp_deriv_l, 4))
})

test_that("gradient of gr_GP_mod_common_hp_k() works for the Linear kernel", {
  k <- seq_len(3)
  m_k <- list("K1" = rep(0, 5), "K2" = rep(0, 5), "K3" = rep(0, 5))

  # db <- simu_db(N = 10, common_input = TRUE)
  # hp_k <- MagmaClustR:::hp("LIN", list_ID = names(m_k))
  # hp_i <- MagmaClustR:::hp("LIN", list_ID = unique(db$ID))


  db <- tibble::tibble("ID" = 1:5, Input = 1:5, Output = 2:6, Covariates = 3:7)
  hp_k <- tibble::tibble("ID" = names(m_k), "lin_slope" = 1, "lin_offset" = 1)
  hp_i <- tibble::tibble("ID" = unique(db$ID), "lin_slope" = 1, "lin_offset" = 1)

  old_hp_mixture <- ini_mixture(data = db, k = length(k), nstart = 50)
  prop_mixture_1 <- old_hp_mixture %>% dplyr::select(-.data$ID)
  hp_k[["prop_mixture"]] <- sapply(
    prop_mixture_1,
    function(x) {
      x %>%
        unlist() %>%
        mean()
    }
  )

  mu_k_param <- MagmaClustR:::ve_step(
    db, m_k, "LIN", "LIN",
    hp_k, hp_i, old_hp_mixture, 2, 0.001
  )

  db_1 <- mu_k_param$mean
  new_cov <- mu_k_param$cov
  hp_i_1 <- hp_i %>%
    dplyr::select(-.data$ID) %>%
    dplyr::slice(1)

  hp_v <- hp_l <- hp_i_1
  hp_v$lin_slope <- hp_i_1$lin_slope + 10^(-8)
  hp_l$lin_offset <- hp_i_1$lin_offset + 10^(-8)

  deriv_v <- gr_GP_mod_common_hp_k(
    hp_i_1, db_1,
    m_k, "LIN", new_cov, 1
  )[["lin_slope"]] %>%
    sum()
  deriv_l <- gr_GP_mod_common_hp_k(
    hp_i_1, db_1, m_k, "LIN",
    new_cov, 1
  )[["lin_offset"]] %>% sum()

  emp_deriv_v <- (elbo_GP_mod_common_hp_k(hp_v, db_1, m_k, "LIN", new_cov, 1) -
    elbo_GP_mod_common_hp_k(hp_i_1, db_1, m_k, "LIN", new_cov, 1)) / 10^(-8)
  emp_deriv_l <- (elbo_GP_mod_common_hp_k(hp_l, db_1, m_k, "LIN", new_cov, 1) -
    elbo_GP_mod_common_hp_k(hp_i_1, db_1, m_k, "LIN", new_cov, 1)) / 10^(-8)

  round(deriv_v, 3) %>% expect_equal(round(as.double(emp_deriv_v), 3))
  round(deriv_l, 3) %>% expect_equal(round(as.double(emp_deriv_l), 3))
})

test_that("gradient of gr_GP_mod_common_hp_k()
          works for the Periodic  kernel", {
  k <- seq_len(3)
  m_k <- list("K1" = rep(0, 5), "K2" = rep(0, 5), "K3" = rep(0, 5))

  # db <- simu_db(N = 10, common_input = TRUE)
  # hp_k <- MagmaClustR:::hp("PERIO", list_ID = names(m_k))
  # hp_i <- MagmaClustR:::hp("PERIO", list_ID = unique(db$ID))


  db <- tibble::tibble("ID" = 1:5, Input = 1:5, Output = 2:6, Covariates = 3:7)
  hp_k <- tibble::tibble(
    "ID" = names(m_k),
    "perio_variance" = 1,
    "perio_lengthscale" = 1,
    "period" = 1
  )
  hp_i <- tibble::tibble(
    "ID" = unique(db$ID),
    "perio_variance" = 1,
    "perio_lengthscale" = 1,
    "period" = 1
  )

  old_hp_mixture <- ini_mixture(data = db, k = length(k), nstart = 50)
  prop_mixture_1 <- old_hp_mixture %>% dplyr::select(-.data$ID)
  hp_k[["prop_mixture"]] <- sapply(
    prop_mixture_1,
    function(x) {
      x %>%
        unlist() %>%
        mean()
    }
  )

  mu_k_param <- MagmaClustR:::ve_step(
    db, m_k, "PERIO", "PERIO",
    hp_k, hp_i, old_hp_mixture, 2, 0.001
  )

  db_1 <- mu_k_param$mean
  new_cov <- mu_k_param$cov
  hp_i_1 <- hp_i %>%
    dplyr::select(-.data$ID) %>%
    dplyr::slice(1)

  hp_v <- hp_l <- hp_p <- hp_i_1
  hp_v$perio_variance <- hp_i_1$perio_variance + 10^(-8)
  hp_l$perio_lengthscale <- hp_i_1$perio_lengthscale + 10^(-8)
  hp_p$period <- hp_i_1$period + 10^(-8)

  deriv_v <- gr_GP_mod_common_hp_k(
    hp_i_1, db_1, m_k, "PERIO",
    new_cov, 1
  )[["perio_variance"]] %>% sum()
  deriv_l <- gr_GP_mod_common_hp_k(
    hp_i_1, db_1, m_k, "PERIO",
    new_cov, 1
  )[["perio_lengthscale"]] %>% sum()
  deriv_p <- gr_GP_mod_common_hp_k(
    hp_i_1, db_1, m_k, "PERIO",
    new_cov, 1
  )[["period"]] %>% sum()

  emp_deriv_v <- (elbo_GP_mod_common_hp_k(hp_v, db_1, m_k, "PERIO", new_cov, 1) -
    elbo_GP_mod_common_hp_k(hp_i_1, db_1, m_k, "PERIO", new_cov, 1)) / 10^(-8)
  emp_deriv_l <- (elbo_GP_mod_common_hp_k(hp_l, db_1, m_k, "PERIO", new_cov, 1) -
    elbo_GP_mod_common_hp_k(hp_i_1, db_1, m_k, "PERIO", new_cov, 1)) / 10^(-8)
  emp_deriv_p <- (elbo_GP_mod_common_hp_k(hp_p, db_1, m_k, "PERIO", new_cov, 1) -
    elbo_GP_mod_common_hp_k(hp_i_1, db_1, m_k, "PERIO", new_cov, 1)) / 10^(-8)

  round(deriv_v, 3) %>% expect_equal(round(as.double(emp_deriv_v), 3))
  round(deriv_l, 2) %>% expect_equal(round(as.double(emp_deriv_l), 2))
  round(deriv_p, 3) %>% expect_equal(round(as.double(emp_deriv_p), 3))
})

test_that("gradient of gr_GP_mod_common_hp_k()
          works for the  Rational Quadratic kernel", {
  k <- seq_len(3)
  m_k <- list("K1" = rep(0, 5), "K2" = rep(0, 5), "K3" = rep(0, 5))

  # db <- simu_db(N = 10, common_input = TRUE)
  # hp_k <- MagmaClustR:::hp("RQ", list_ID = names(m_k))
  # hp_i <- MagmaClustR:::hp("RQ", list_ID = unique(db$ID))


  db <- tibble::tibble("ID" = 1:5, Input = 1:5, Output = 2:6, Covariates = 3:7)
  hp_k <- tibble::tibble(
    "ID" = names(m_k),
    "rq_variance" = 1,
    "rq_lengthscale" = 1,
    "rq_scale" = 1
  )
  hp_i <- tibble::tibble(
    "ID" = unique(db$ID),
    "rq_variance" = 1,
    "rq_lengthscale" = 1,
    "rq_scale" = 1
  )

  old_hp_mixture <- ini_mixture(data = db, k = length(k), nstart = 50)
  prop_mixture_1 <- old_hp_mixture %>% dplyr::select(-.data$ID)
  hp_k[["prop_mixture"]] <- sapply(
    prop_mixture_1,
    function(x) {
      x %>%
        unlist() %>%
        mean()
    }
  )

  mu_k_param <- MagmaClustR:::ve_step(
    db, m_k, "RQ", "RQ",
    hp_k, hp_i, old_hp_mixture,2, 0.001
  )

  db_1 <- mu_k_param$mean
  new_cov <- mu_k_param$cov
  hp_i_1 <- hp_i %>%
    dplyr::select(-.data$ID) %>%
    dplyr::slice(1)

  hp_v <- hp_l <- hp_p <- hp_i_1
  hp_v$rq_variance <- hp_i_1$rq_variance + 10^(-8)
  hp_l$rq_lengthscale <- hp_i_1$rq_lengthscale + 10^(-8)
  hp_p$rq_scale <- hp_i_1$rq_scale + 10^(-8)

  deriv_v <- gr_GP_mod_common_hp_k(
    hp_i_1, db_1, m_k, "RQ",
    new_cov, 1
  )[["rq_variance"]] %>% sum()
  deriv_l <- gr_GP_mod_common_hp_k(
    hp_i_1, db_1, m_k, "RQ",
    new_cov, 1
  )[["rq_lengthscale"]] %>% sum()
  deriv_p <- gr_GP_mod_common_hp_k(
    hp_i_1, db_1, m_k, "RQ",
    new_cov, 1
  )[["rq_scale"]] %>% sum()

  emp_deriv_v <- (elbo_GP_mod_common_hp_k(hp_v, db_1, m_k, "RQ", new_cov, 1) -
    elbo_GP_mod_common_hp_k(hp_i_1, db_1, m_k, "RQ", new_cov, 1)) / 10^(-8)
  emp_deriv_l <- (elbo_GP_mod_common_hp_k(hp_l, db_1, m_k, "RQ", new_cov, 1) -
    elbo_GP_mod_common_hp_k(hp_i_1, db_1, m_k, "RQ", new_cov, 1)) / 10^(-8)
  emp_deriv_p <- (elbo_GP_mod_common_hp_k(hp_p, db_1, m_k, "RQ", new_cov, 1) -
    elbo_GP_mod_common_hp_k(hp_i_1, db_1, m_k, "RQ", new_cov, 1)) / 10^(-8)

  round(deriv_v, 3) %>% expect_equal(round(as.double(emp_deriv_v), 3))
  round(deriv_l, 3) %>% expect_equal(round(as.double(emp_deriv_l), 3))
  round(deriv_p, 3) %>% expect_equal(round(as.double(emp_deriv_p), 3))
})

########################## gr_clust_multi_GP_common_hp_i() #####################

test_that("gradient of gr_clust_multi_GP_common_hp_i()
          works for the Squared Exponential kernel", {
  k <- seq_len(3)
  m_k <- list("K1" = rep(0, 5), "K2" = rep(0, 5), "K3" = rep(0, 5))

  # db <- simu_db(N = 10, common_input = TRUE)
  # hp_k <- MagmaClustR:::hp("SE", list_ID = names(m_k))
  # hp_i <- MagmaClustR:::hp("SE", list_ID = unique(db$ID))


  db <- tibble::tibble("ID" = 1:5, Input = 1:5, Output = 2:6, Covariates = 3:7)
  hp_k <- tibble::tibble(
    "ID" = names(m_k),
    "se_variance" = 1,
    "se_lengthscale" = 1
  )
  hp_i <- tibble::tibble(
    "ID" = unique(db$ID),
    "se_variance" = 1,
    "se_lengthscale" = 1
  )

  old_hp_mixture <- ini_mixture(data = db, k = length(k), nstart = 50)
  prop_mixture_1 <- old_hp_mixture %>% dplyr::select(-.data$ID)
  hp_k[["prop_mixture"]] <- sapply(
    prop_mixture_1,
    function(x) {
      x %>%
        unlist() %>%
        mean()
    }
  )

  mu_k_param <- MagmaClustR:::ve_step(
    db, m_k, "SE", "SE", hp_k, hp_i,
    old_hp_mixture, 2, 0.001
  )

  db_1 <- db %>% dplyr::filter(.data$ID == 1)
  hp_i_1 <- hp_i %>%
    dplyr::select(-.data$ID) %>%
    dplyr::slice(1)

  hp_v <- hp_l <- hp_i_1
  hp_v$se_variance <- hp_i_1$se_variance + 10^(-8)
  hp_l$se_lengthscale <- hp_i_1$se_lengthscale + 10^(-8)



  deriv_v <- gr_clust_multi_GP_common_hp_i(
    hp_i_1, db_1, mu_k_param,
    "SE", 1
  )[["se_variance"]]
  deriv_l <- gr_clust_multi_GP_common_hp_i(
    hp_i_1, db_1, mu_k_param,
    "SE", 1
  )[["se_lengthscale"]]

  emp_deriv_v <- (elbo_clust_multi_GP_common_hp_i(
    hp_v, db_1,
    mu_k_param, "SE", 1
  ) -
    elbo_clust_multi_GP_common_hp_i(hp_i_1, db_1, mu_k_param, "SE", 1)) / 10^(-8)
  emp_deriv_l <- (elbo_clust_multi_GP_common_hp_i(
    hp_l, db_1,
    mu_k_param, "SE", 1
  ) -
    elbo_clust_multi_GP_common_hp_i(hp_i_1, db_1, mu_k_param, "SE", 1)) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(as.double(emp_deriv_v), 4))
  round(deriv_l, 4) %>% expect_equal(round(as.double(emp_deriv_l), 4))
})

test_that("gradient of gr_clust_multi_GP_common_hp_i()
          works for the Linear kernel", {
  k <- seq_len(3)
  m_k <- list("K1" = rep(0, 5), "K2" = rep(0, 5), "K3" = rep(0, 5))

  # db <- simu_db(N = 10, common_input = TRUE)
  # hp_k <- MagmaClustR:::hp("LIN", list_ID = names(m_k))
  # hp_i <- MagmaClustR:::hp("LIN", list_ID = unique(db$ID))


  db <- tibble::tibble("ID" = 1:5, Input = 1:5, Output = 2:6, Covariates = 3:7)
  hp_k <- tibble::tibble(
    "ID" = names(m_k),
    "lin_slope" = 1, "lin_offset" = 1
  )
  hp_i <- tibble::tibble(
    "ID" = unique(db$ID),
    "lin_slope" = 1, "lin_offset" = 1
  )

  old_hp_mixture <- ini_mixture(data = db, k = length(k), nstart = 50)
  prop_mixture_1 <- old_hp_mixture %>% dplyr::select(-.data$ID)
  hp_k[["prop_mixture"]] <- sapply(
    prop_mixture_1,
    function(x) {
      x %>%
        unlist() %>%
        mean()
    }
  )

  mu_k_param <- MagmaClustR:::ve_step(
    db, m_k, "LIN", "LIN",
    hp_k, hp_i, old_hp_mixture, 2, 0.001
  )

  db_1 <- db %>% dplyr::filter(.data$ID == 1)
  hp_i_1 <- hp_i %>%
    dplyr::select(-.data$ID) %>%
    dplyr::slice(1)

  hp_v <- hp_l <- hp_i_1
  hp_v$lin_slope <- hp_i_1$lin_slope + 10^(-8)
  hp_l$lin_offset <- hp_i_1$lin_offset + 10^(-8)



  deriv_v <- gr_clust_multi_GP_common_hp_i(
    hp_i_1, db_1, mu_k_param,
    "LIN", 1
  )[["lin_slope"]]
  deriv_l <- gr_clust_multi_GP_common_hp_i(
    hp_i_1, db_1, mu_k_param,
    "LIN", 1
  )[["lin_offset"]]

  emp_deriv_v <- (elbo_clust_multi_GP_common_hp_i(
    hp_v, db_1,
    mu_k_param, "LIN", 1
  ) -
    elbo_clust_multi_GP_common_hp_i(hp_i_1, db_1, mu_k_param, "LIN", 1)) / 10^(-8)
  emp_deriv_l <- (elbo_clust_multi_GP_common_hp_i(
    hp_l, db_1,
    mu_k_param, "LIN", 1
  ) -
    elbo_clust_multi_GP_common_hp_i(hp_i_1, db_1, mu_k_param, "LIN", 1)) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(as.double(emp_deriv_v), 4))
  round(deriv_l, 4) %>% expect_equal(round(as.double(emp_deriv_l), 4))
})

test_that("gradient of gr_clust_multi_GP_common_hp_i()
          works for the Linear kernel", {
  k <- seq_len(3)
  m_k <- list("K1" = rep(0, 5), "K2" = rep(0, 5), "K3" = rep(0, 5))

  # db <- simu_db(N = 10, common_input = TRUE)
  # hp_k <- MagmaClustR:::hp("PERIO", list_ID = names(m_k))
  # hp_i <- MagmaClustR:::hp("PERIO", list_ID = unique(db$ID))


  db <- tibble::tibble("ID" = 1:5, Input = 1:5, Output = 2:6, Covariates = 3:7)
  hp_k <- tibble::tibble(
    "ID" = names(m_k),
    "perio_variance" = 1,
    "perio_lengthscale" = 1,
    "period" = 1
  )
  hp_i <- tibble::tibble(
    "ID" = unique(db$ID),
    "perio_variance" = 1,
    "perio_lengthscale" = 1,
    "period" = 1
  )

  old_hp_mixture <- ini_mixture(data = db, k = length(k), nstart = 50)
  prop_mixture_1 <- old_hp_mixture %>% dplyr::select(-.data$ID)
  hp_k[["prop_mixture"]] <- sapply(
    prop_mixture_1,
    function(x) {
      x %>%
        unlist() %>%
        mean()
    }
  )

  mu_k_param <- MagmaClustR:::ve_step(
    db, m_k, "PERIO", "PERIO",
    hp_k, hp_i, old_hp_mixture, 2, 0.001
  )

  db_1 <- db %>% dplyr::filter(.data$ID == 1)
  hp_i_1 <- hp_i %>%
    dplyr::select(-.data$ID) %>%
    dplyr::slice(1)

  hp_v <- hp_l <- hp_p <- hp_i_1
  hp_v$perio_variance <- hp_i_1$perio_variance + 10^(-8)
  hp_l$perio_lengthscale <- hp_i_1$perio_lengthscale + 10^(-8)
  hp_p$period <- hp_i_1$period + 10^(-8)


  deriv_v <- gr_clust_multi_GP_common_hp_i(
    hp_i_1, db_1, mu_k_param,
    "PERIO", 1
  )[["perio_variance"]]
  deriv_l <- gr_clust_multi_GP_common_hp_i(
    hp_i_1, db_1, mu_k_param,
    "PERIO", 1
  )[["perio_lengthscale"]]
  deriv_p <- gr_clust_multi_GP_common_hp_i(
    hp_i_1, db_1, mu_k_param,
    "PERIO", 1
  )[["period"]]


  emp_deriv_v <- (elbo_clust_multi_GP_common_hp_i(
    hp_v, db_1,
    mu_k_param, "PERIO", 1
  ) -
    elbo_clust_multi_GP_common_hp_i(hp_i_1, db_1, mu_k_param, "PERIO", 1)) / 10^(-8)
  emp_deriv_l <- (elbo_clust_multi_GP_common_hp_i(
    hp_l, db_1,
    mu_k_param, "PERIO", 1
  ) -
    elbo_clust_multi_GP_common_hp_i(hp_i_1, db_1, mu_k_param, "PERIO", 1)) / 10^(-8)
  emp_deriv_p <- (elbo_clust_multi_GP_common_hp_i(
    hp_p, db_1,
    mu_k_param, "PERIO", 1
  ) -
    elbo_clust_multi_GP_common_hp_i(hp_i_1, db_1, mu_k_param, "PERIO", 1)) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(as.double(emp_deriv_v), 4))
  round(deriv_l, 4) %>% expect_equal(round(as.double(emp_deriv_l), 4))
})

test_that("gradient of gr_clust_multi_GP_common_hp_i()
          works for the Rational Quadratic kernel", {
  k <- seq_len(3)
  m_k <- list("K1" = rep(0, 5), "K2" = rep(0, 5), "K3" = rep(0, 5))

  # db <- simu_db(N = 10, common_input = TRUE)
  # hp_k <- MagmaClustR:::hp("RQ", list_ID = names(m_k))
  # hp_i <- MagmaClustR:::hp("RQ", list_ID = unique(db$ID))


  db <- tibble::tibble("ID" = 1:5, Input = 1:5, Output = 2:6, Covariates = 3:7)
  hp_k <- tibble::tibble(
    "ID" = names(m_k),
    "rq_variance" = 1,
    "rq_lengthscale" = 1,
    "rq_scale" = 1
  )
  hp_i <- tibble::tibble(
    "ID" = unique(db$ID),
    "rq_variance" = 1,
    "rq_lengthscale" = 1,
    "rq_scale" = 1
  )


  old_hp_mixture <- ini_mixture(data = db, k = length(k), nstart = 50)
  prop_mixture_1 <- old_hp_mixture %>% dplyr::select(-.data$ID)
  hp_k[["prop_mixture"]] <- sapply(
    prop_mixture_1,
    function(x) {
      x %>%
        unlist() %>%
        mean()
    }
  )

  mu_k_param <- MagmaClustR:::ve_step(
    db, m_k, "RQ", "RQ", hp_k,
    hp_i, old_hp_mixture, 2, 0.001
  )

  db_1 <- db %>% dplyr::filter(.data$ID == 1)
  hp_i_1 <- hp_i %>%
    dplyr::select(-.data$ID) %>%
    dplyr::slice(1)

  hp_v <- hp_l <- hp_p <- hp_i_1
  hp_v$rq_variance <- hp_i_1$rq_variance + 10^(-8)
  hp_l$rq_lengthscale <- hp_i_1$rq_lengthscale + 10^(-8)
  hp_p$rq_scale <- hp_i_1$rq_scale + 10^(-8)


  deriv_v <- gr_clust_multi_GP_common_hp_i(
    hp_i_1, db_1, mu_k_param,
    "RQ", 1
  )[["rq_variance"]]
  deriv_l <- gr_clust_multi_GP_common_hp_i(
    hp_i_1, db_1, mu_k_param,
    "RQ", 1
  )[["rq_lengthscale"]]
  deriv_p <- gr_clust_multi_GP_common_hp_i(
    hp_i_1, db_1, mu_k_param,
    "RQ", 1
  )[["rq_scale"]]


  emp_deriv_v <- (elbo_clust_multi_GP_common_hp_i(
    hp_v, db_1,
    mu_k_param, "RQ", 1
  ) -
    elbo_clust_multi_GP_common_hp_i(hp_i_1, db_1, mu_k_param, "RQ", 1)) / 10^(-8)
  emp_deriv_l <- (elbo_clust_multi_GP_common_hp_i(
    hp_l, db_1,
    mu_k_param, "RQ", 1
  ) -
    elbo_clust_multi_GP_common_hp_i(hp_i_1, db_1, mu_k_param, "RQ", 1)) / 10^(-8)
  emp_deriv_p <- (elbo_clust_multi_GP_common_hp_i(
    hp_p, db_1,
    mu_k_param, "RQ", 1
  ) -
    elbo_clust_multi_GP_common_hp_i(hp_i_1, db_1, mu_k_param, "RQ", 1)) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(as.double(emp_deriv_v), 4))
  round(deriv_l, 4) %>% expect_equal(round(as.double(emp_deriv_l), 4))
})
