test_that("gradient of gr_GP() works for Squared Exponential kernel", {
  db <- tibble::tibble(
    Input = 1:5,
    Output = 2:6,
    Covariate = 3:7,
    Reference = paste(Input, Covariate, sep = ':')
  )
  mean <- rep(0, 5)
  hp <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5)
  new_cov <- kern_to_cov(db %>% dplyr::select(-Output), "SE", hp)

  hp_v <- tibble::tibble(se_variance = 1 + 10^(-8), se_lengthscale = 0.5)
  hp_l <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5 + 10^(-8))

  deriv_v <- gr_GP(hp, db, mean, "SE", new_cov, 0)[["se_variance"]]
  deriv_l <- gr_GP(hp, db, mean, "SE", new_cov, 0)[["se_lengthscale"]]

  emp_deriv_v <- (logL_GP(hp_v, db, mean, "SE", new_cov, 0) -
    logL_GP(hp, db, mean, "SE", new_cov, 0)) /
    10^(-8)
  emp_deriv_l <- (logL_GP(hp_l, db, mean, "SE", new_cov, 0) -
    logL_GP(hp, db, mean, "SE", new_cov, 0)) /
    10^(-8)

  round(deriv_v, 3) %>% expect_equal(round(emp_deriv_v, 3))
  round(deriv_l, 3) %>% expect_equal(round(emp_deriv_l, 3))
})

## TODO: test for gr_GP() for the RQ and PERIOD kernels

test_that("gradient of logL_GP_mod() works for Squared Exponential kernel", {
  db <- tibble::tibble(
    Input = 1:5,
    Output = 2:6,
    Covariate = 3:7,
    Reference = paste(Input, Covariate, sep = ':')
  )

  mean <- rep(0, 5)
  hp <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5)
  new_cov <- kern_to_cov(db %>% dplyr::select(-Output), "SE", hp)

  hp_v <- tibble::tibble(se_variance = 1 + 10^(-8), se_lengthscale = 0.5)
  hp_l <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5 + 10^(-8))

  deriv_v <- gr_GP_mod(hp, db, mean, "SE", new_cov, 0)[["se_variance"]]
  deriv_l <- gr_GP_mod(hp, db, mean, "SE", new_cov, 0)[["se_lengthscale"]]

  emp_deriv_v <- (logL_GP_mod(hp_v, db, mean, "SE", new_cov, 0) -
    logL_GP_mod(hp, db, mean, "SE", new_cov, 0)) /
    10^(-8)
  emp_deriv_l <- (logL_GP_mod(hp_l, db, mean, "SE", new_cov, 0) -
    logL_GP_mod(hp, db, mean, "SE", new_cov, 0)) /
    10^(-8)

  round(deriv_v, 3) %>% expect_equal(round(emp_deriv_v, 3))
  round(deriv_l, 3) %>% expect_equal(round(emp_deriv_l, 3))
})

## TODO: test for gr_GP_mod() for the RQ and PERIOD kernels

test_that("gradient of logL_GP_mod_shared_hp() works", {
  db <- tibble::tibble(
    ID = rep(1:5, each = 4),
    Output = 1:20,
    Input = 2:21,
    Covariate = c(1:10, 23, 77, 1:8),
    Reference = paste(Input, Covariate, sep = ':')
  )
  mean <- tibble::tibble(
    "Input" = db$Input,
    "Covariate" = db$Covariate,
    "Reference" = db$Reference,
    "Output" = 0
  )
  hp <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5)
  new_cov <- kern_to_cov(db %>% dplyr::select(-Output), "SE", hp)

  hp_v <- tibble::tibble(se_variance = 1 + 10^(-8), se_lengthscale = 0.5)
  hp_l <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5 + 10^(-8))

  deriv_v <- gr_GP_mod_shared_hp(
    hp,
    db,
    mean,
    "SE",
    new_cov,
    0.1
  )[["se_variance"]]
  deriv_l <- gr_GP_mod_shared_hp(
    hp,
    db,
    mean,
    "SE",
    new_cov,
    0
  )[["se_lengthscale"]]

  emp_deriv_v <- (logL_GP_mod_shared_hp(hp_v, db, mean, "SE", new_cov, 0.1) -
    logL_GP_mod_shared_hp(hp, db, mean, "SE", new_cov, 0.1)) /
    10^(-8)

  emp_deriv_l <- (logL_GP_mod_shared_hp(hp_l, db, mean, "SE", new_cov, 0) -
    logL_GP_mod_shared_hp(hp, db, mean, "SE", new_cov, 0)) /
    10^(-8)

  round(deriv_v, 3) %>% expect_equal(round(emp_deriv_v, 3))
  round(deriv_l, 3) %>% expect_equal(round(emp_deriv_l, 3))
})

## TODO: test for gr_GP_mod_shared_hp() for the RQ and PERIOD kernels

test_that("gradient of gr_sum_logL_GP_clust() works for SE kernel", {
  db <- tibble::tibble(
    Input = 1:5,
    Output = 2:6,
    Covariate = 3:7,
    Reference = paste(1:5, 3:7, sep = ':')
  )
  mean <- list(
    "K1" = tibble::tibble(
      "Input" = db$Input,
      "Covariate" = db$Covariate,
      "Reference" = db$Reference,
      "Output" = 0
    ),
    "K2" = tibble::tibble(
      "Input" = db$Input,
      "Covariate" = db$Covariate,
      "Reference" = db$Reference,
      "Output" = 0
    )
  )
  hp <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5)
  new_cov <- list(
    "K1" = kern_to_cov(db %>% dplyr::select(-Output), "SE", hp),
    "K2" = kern_to_cov(db %>% dplyr::select(-Output), "SE", hp)
  )
  mixture <- tibble::tibble("K1" = 0.4, "K2" = 0.6)
  hp_v <- tibble::tibble(se_variance = 1 + 10^(-8), se_lengthscale = 0.5)
  hp_l <- tibble::tibble(se_variance = 1, se_lengthscale = 0.5 + 10^(-8))

  deriv_v <- gr_sum_logL_GP_clust(
    hp,
    db,
    mixture,
    mean,
    "SE",
    new_cov,
    0
  )[["se_variance"]]
  deriv_l <- gr_sum_logL_GP_clust(
    hp,
    db,
    mixture,
    mean,
    "SE",
    new_cov,
    0
  )[["se_lengthscale"]]

  emp_deriv_v <- (sum_logL_GP_clust(
    hp_v,
    db,
    mixture,
    mean,
    "SE",
    new_cov,
    NULL,
    0
  ) -
    sum_logL_GP_clust(hp, db, mixture, mean, "SE", new_cov, NULL, 0)) /
    10^(-8)
  emp_deriv_l <- (sum_logL_GP_clust(
    hp_l,
    db,
    mixture,
    mean,
    "SE",
    new_cov,
    NULL,
    0
  ) -
    sum_logL_GP_clust(hp, db, mixture, mean, "SE", new_cov, NULL, 0)) /
    10^(-8)

  round(deriv_v, 3) %>% expect_equal(round(emp_deriv_v, 3))
  round(deriv_l, 3) %>% expect_equal(round(emp_deriv_l, 3))
})

## ---- Multi-output convolution kernel gradients --------------------------

test_that("gradient of gr_GP_mod() works for the multi-output convolution kernel", {
  set.seed(6274)
  inputs <- tibble::tibble(
    Output_ID = factor(rep(1:2, each = 3)),
    Input = rep(c(0.2, 0.5, 0.8), 2)
  )
  inputs$Reference <- paste0("o", inputs$Output_ID, ";", inputs$Input)
  n <- nrow(inputs)
  db <- inputs
  db$Output <- stats::rnorm(n)
  mean_vec <- rep(0, n)
  A <- matrix(stats::rnorm(n * n), n)
  post_cov <- crossprod(A) / n * 0.05

  hp_tbl <- hp(convolution_kernel, list_task_ID = "a", list_output_ID = 1:2,
               shared_hp_outputs = FALSE, noise = TRUE) %>%
    dplyr::select(-Task_ID)
  flat <- flatten_hp_mo(hp_tbl)

  ana <- gr_GP_mod(flat, db, mean_vec, convolution_kernel, post_cov, 1e-10)
  fd <- sapply(names(flat), function(p) {
    eps <- 1e-6
    hp_p <- flat; hp_m <- flat
    hp_p[p] <- hp_p[p] + eps; hp_m[p] <- hp_m[p] - eps
    (logL_GP_mod(hp_p, db, mean_vec, convolution_kernel, post_cov, 1e-10) -
     logL_GP_mod(hp_m, db, mean_vec, convolution_kernel, post_cov, 1e-10)) / (2 * eps)
  })
  expect_equal(as.numeric(ana), as.numeric(fd), tolerance = 1e-4)
})

test_that("gradient of gr_GP_mod_shared_hp() works for the multi-output kernel", {
  set.seed(7731)
  grid <- tibble::tibble(Output_ID = factor(rep(1:2, each = 3)),
                         Input = rep(c(0.2, 0.5, 0.8), 2))
  grid$Reference <- paste0("o", grid$Output_ID, "_", grid$Input)
  grid <- dplyr::arrange(grid, .data$Reference)
  refs <- grid$Reference
  db <- dplyr::bind_rows(
    dplyr::mutate(grid, ID = "i1", Output = stats::rnorm(6)),
    dplyr::mutate(grid, ID = "i2", Output = stats::rnorm(6))
  )
  db <- dplyr::arrange(db, .data$ID, .data$Reference)
  mean_tbl <- tibble::tibble(Reference = refs, Output = stats::rnorm(6))
  A <- matrix(stats::rnorm(36), 6)
  post_cov <- crossprod(A) / 6 * 0.05
  rownames(post_cov) <- colnames(post_cov) <- refs

  hp_tbl <- hp(convolution_kernel, list_task_ID = "a", list_output_ID = 1:2,
               shared_hp_outputs = FALSE, noise = TRUE) %>%
    dplyr::select(-Task_ID)
  flat <- flatten_hp_mo(hp_tbl)

  ana <- gr_GP_mod_shared_hp(flat, db, mean_tbl, convolution_kernel, post_cov, 1e-10)
  fd <- sapply(names(flat), function(p) {
    eps <- 1e-6
    hp_p <- flat; hp_m <- flat
    hp_p[p] <- hp_p[p] + eps; hp_m[p] <- hp_m[p] - eps
    (logL_GP_mod_shared_hp(hp_p, db, mean_tbl, convolution_kernel, post_cov, 1e-10) -
     logL_GP_mod_shared_hp(hp_m, db, mean_tbl, convolution_kernel, post_cov, 1e-10)) / (2 * eps)
  })
  expect_equal(as.numeric(ana), as.numeric(fd), tolerance = 1e-4)
})
