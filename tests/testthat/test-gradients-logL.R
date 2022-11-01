test_that("gradient of gr_GP() works for Squared Exponential kernel", {
  db <- tibble::tibble(Input = 1:5, Output = 2:6, Covariate = 3:7,
                       Reference = paste(Input, Covariate, sep = ':'))
  mean <- rep(0, 5)
  hp <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5)
  new_cov <- kern_to_cov(db %>% dplyr::select(- Output), "SE", hp)

  hp_v <- tibble::tibble(se_variance = 1 + 10^(-8), se_lengthscale = 0.5)
  hp_l <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5 + 10^(-8))

  deriv_v <- gr_GP(hp, db, mean, "SE", new_cov, 0)[["se_variance"]]
  deriv_l <- gr_GP(hp, db, mean, "SE", new_cov, 0)[["se_lengthscale"]]

  emp_deriv_v <- (logL_GP(hp_v, db, mean, "SE", new_cov, 0) -
    logL_GP(hp, db, mean, "SE", new_cov, 0)) / 10^(-8)
  emp_deriv_l <- (logL_GP(hp_l, db, mean, "SE", new_cov, 0) -
    logL_GP(hp, db, mean, "SE", new_cov, 0)) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(emp_deriv_v, 4))
  round(deriv_l, 4) %>% expect_equal(round(emp_deriv_l, 4))
})

## TODO: test for gr_GP() for the RQ and PERIOD kernels

test_that("gradient of logL_GP_mod() works for Squared Exponential kernel", {
  db <- tibble::tibble(Input = 1:5, Output = 2:6, Covariate = 3:7,
                       Reference = paste(Input, Covariate, sep = ':'))

  mean <- rep(0, 5)
  hp <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5)
  new_cov <- kern_to_cov(db %>% dplyr::select(- Output), "SE", hp)

  hp_v <- tibble::tibble(se_variance = 1 + 10^(-8), se_lengthscale = 0.5)
  hp_l <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5 + 10^(-8))

  deriv_v <- gr_GP_mod(hp, db, mean, "SE", new_cov, 0)[["se_variance"]]
  deriv_l <- gr_GP_mod(hp, db, mean, "SE", new_cov, 0)[["se_lengthscale"]]

  emp_deriv_v <- (logL_GP_mod(hp_v, db, mean, "SE", new_cov, 0) -
    logL_GP_mod(hp, db, mean, "SE", new_cov, 0)) / 10^(-8)
  emp_deriv_l <- (logL_GP_mod(hp_l, db, mean, "SE", new_cov, 0) -
    logL_GP_mod(hp, db, mean, "SE", new_cov, 0)) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(emp_deriv_v, 4))
  round(deriv_l, 4) %>% expect_equal(round(emp_deriv_l, 4))
})

## TODO: test for gr_GP_mod() for the RQ and PERIOD kernels

test_that("gradient of logL_GP_mod_common_hp() works", {
  db <- tibble::tibble(
    ID = rep(1:5, each = 4),
    Output = 1:20,
    Input = 2:21,
    Covariate = c(1:10, 23, 77, 1:8),
    Reference = paste(Input, Covariate, sep = ':')
  )
  mean <- tibble::tibble("Input" = db$Input, "Covariate" = db$Covariate,
                         "Reference" = db$Reference, "Output" = 0)
  hp <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5)
  new_cov <- kern_to_cov(db %>% dplyr::select(- Output), "SE", hp)

  hp_v <- tibble::tibble(se_variance = 1 + 10^(-8), se_lengthscale = 0.5)
  hp_l <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5 + 10^(-8))

  deriv_v <- gr_GP_mod_common_hp(
    hp, db, mean,
    "SE", new_cov, 0.1
  )[["se_variance"]]
  deriv_l <- gr_GP_mod_common_hp(
    hp, db, mean,
    "SE", new_cov, 0
  )[["se_lengthscale"]]

  emp_deriv_v <- (logL_GP_mod_common_hp(hp_v, db, mean, "SE", new_cov, 0.1) -
    logL_GP_mod_common_hp(hp, db, mean, "SE", new_cov, 0.1)) / 10^(-8)

  emp_deriv_l <- (logL_GP_mod_common_hp(hp_l, db, mean, "SE", new_cov, 0) -
    logL_GP_mod_common_hp(hp, db, mean, "SE", new_cov, 0)) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(emp_deriv_v, 4))
  round(deriv_l, 4) %>% expect_equal(round(emp_deriv_l, 4))
})

## TODO: test for gr_GP_mod_common_hp() for the RQ and PERIOD kernels

test_that("gradient of gr_sum_logL_GP_clust() works for SE kernel", {
  db <- tibble::tibble(Input = 1:5, Output = 2:6,
                       Covariate = 3:7, Reference = paste(1:5, 3:7, sep = ':'))
  mean <- list(
    "K1" = tibble::tibble("Input" = db$Input, "Covariate" = db$Covariate,
                           "Reference" = db$Reference, "Output" = 0),
    "K2" = tibble::tibble("Input" = db$Input, "Covariate" = db$Covariate,
                                "Reference" = db$Reference, "Output" = 0)
  )
  hp <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5)
  new_cov <- list(
    "K1" = kern_to_cov(db %>% dplyr::select(- Output), "SE", hp),
    "K2" = kern_to_cov(db %>% dplyr::select(- Output), "SE", hp)
  )
  mixture <- tibble::tibble("K1" = 0.4, "K2" = 0.6)
  hp_v <- tibble::tibble(se_variance = 1 + 10^(-8), se_lengthscale = 0.5)
  hp_l <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5 + 10^(-8))

  deriv_v <- gr_sum_logL_GP_clust(
    hp, db, mixture, mean,
    "SE", new_cov, 0
  )[["se_variance"]]
  deriv_l <- gr_sum_logL_GP_clust(
    hp, db, mixture, mean,
    "SE", new_cov, 0
  )[["se_lengthscale"]]

  emp_deriv_v <- (sum_logL_GP_clust(
    hp_v, db, mixture, mean,
    "SE", new_cov, NULL, 0
  ) -
    sum_logL_GP_clust(hp, db, mixture, mean, "SE", new_cov, NULL, 0)) / 10^(-8)
  emp_deriv_l <- (sum_logL_GP_clust(
    hp_l, db, mixture, mean,
    "SE", new_cov, NULL, 0
  ) -
    sum_logL_GP_clust(hp, db, mixture, mean, "SE", new_cov, NULL, 0)) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(emp_deriv_v, 4))
  round(deriv_l, 4) %>% expect_equal(round(emp_deriv_l, 4))
})
