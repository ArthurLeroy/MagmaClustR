test_that("Squared Exponential kernel works for null distance", {
  hp <- tibble::tibble(se_variance = 1, se_lengthscale = 1)

  se_kernel(2, 2, hp)[1] %>% expect_equal(exp(1))
})

test_that("Periodic kernel works for null distance", {
  hp <- tibble::tibble(perio_variance = 1, perio_lengthscale = 1, period = 1)

  perio_kernel(2, 2, hp)[1] %>% expect_equal(exp(1))
})

test_that("Rational quadratic kernel works for null distance", {
  hp <- tibble::tibble(rq_variance = 1, rq_lengthscale = 1, rq_scale = 1)

  rq_kernel(2, 2, hp)[1] %>% expect_equal(exp(1))
})

test_that("Linear kernel works for null distance", {
  hp <- tibble::tibble(lin_slope = 1, lin_offset = 1)

  lin_kernel(0, 2, hp)[1] %>% expect_equal(exp(1))
})

test_that("gradients for the Squared Exponential kernel are valid", {
  hp <- tibble::tibble(se_variance = 1, se_lengthscale = 1)
  hp_v <- tibble::tibble(se_variance = 1 + 10^(-8), se_lengthscale = 1)
  hp_l <- tibble::tibble(se_variance = 1, se_lengthscale = 1 + 10^(-8))

  deriv_v <- se_kernel(c(1, 2), c(2, 3), hp, "se_variance")
  deriv_l <- se_kernel(c(1, 2), c(2, 3), hp, "se_lengthscale")

  emp_deriv_v <- (se_kernel(c(1, 2), c(2, 3), hp_v)[1] -
    se_kernel(c(1, 2), c(2, 3), hp)[1]) /
    10^(-8)
  emp_deriv_l <- (se_kernel(c(1, 2), c(2, 3), hp_l)[1] -
    se_kernel(c(1, 2), c(2, 3), hp)[1]) /
    10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(emp_deriv_v, 4))
  round(deriv_l, 4) %>% expect_equal(round(emp_deriv_l, 4))
})

test_that("gradients for the Periodic kernel are valid", {
  hp <- tibble::tibble(
    perio_variance = 1,
    perio_lengthscale = 1,
    period = pi
  )
  hp_v <- tibble::tibble(
    perio_variance = 1 + 10^(-8),
    perio_lengthscale = 1,
    period = pi
  )
  hp_l <- tibble::tibble(
    perio_variance = 1,
    perio_lengthscale = 1 + 10^(-8),
    period = pi
  )
  hp_p <- tibble::tibble(
    perio_variance = 1,
    perio_lengthscale = 1,
    period = pi + 10^(-8)
  )

  deriv_v <- perio_kernel(c(1, 2), c(2, 3), hp, "perio_variance")
  deriv_l <- perio_kernel(c(1, 2), c(2, 3), hp, "perio_lengthscale")
  deriv_p <- perio_kernel(c(1, 2), c(2, 3), hp, "period")

  emp_deriv_v <- (perio_kernel(c(1, 2), c(2, 3), hp_v)[1] -
    perio_kernel(c(1, 2), c(2, 3), hp)[1]) /
    10^(-8)
  emp_deriv_l <- (perio_kernel(c(1, 2), c(2, 3), hp_l)[1] -
    perio_kernel(c(1, 2), c(2, 3), hp)[1]) /
    10^(-8)
  emp_deriv_p <- (perio_kernel(c(1, 2), c(2, 3), hp_p)[1] -
    perio_kernel(c(1, 2), c(2, 3), hp)[1]) /
    10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(emp_deriv_v, 4))
  round(deriv_l, 4) %>% expect_equal(round(emp_deriv_l, 4))
  round(deriv_p, 4) %>% expect_equal(round(emp_deriv_p, 4))
})

test_that("gradients for the Rational Quadratic kernel are valid", {
  hp <- tibble::tibble(rq_variance = 1, rq_lengthscale = 1, rq_scale = 1)
  hp_v <- tibble::tibble(
    rq_variance = 1 + 10^(-8),
    rq_lengthscale = 1,
    rq_scale = 1
  )
  hp_l <- tibble::tibble(
    rq_variance = 1,
    rq_lengthscale = 1 + 10^(-8),
    rq_scale = 1
  )
  hp_s <- tibble::tibble(
    rq_variance = 1,
    rq_lengthscale = 1,
    rq_scale = 1 + 10^(-8)
  )

  deriv_v <- rq_kernel(c(1, 2), c(2, 3), hp, "rq_variance")
  deriv_l <- rq_kernel(c(1, 2), c(2, 3), hp, "rq_lengthscale")
  deriv_s <- rq_kernel(c(1, 2), c(2, 3), hp, "rq_scale")

  emp_deriv_v <- (rq_kernel(c(1, 2), c(2, 3), hp_v)[1] -
    rq_kernel(c(1, 2), c(2, 3), hp)[1]) /
    10^(-8)
  emp_deriv_l <- (rq_kernel(c(1, 2), c(2, 3), hp_l)[1] -
    rq_kernel(c(1, 2), c(2, 3), hp)[1]) /
    10^(-8)
  emp_deriv_s <- (rq_kernel(c(1, 2), c(2, 3), hp_s)[1] -
    rq_kernel(c(1, 2), c(2, 3), hp)[1]) /
    10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(emp_deriv_v, 4))
  round(deriv_l, 4) %>% expect_equal(round(emp_deriv_l, 4))
  round(deriv_s, 4) %>% expect_equal(round(emp_deriv_s, 4))
})

test_that("gradients for the Linear kernel are valid", {
  hp <- tibble::tibble(lin_slope = 1, lin_intercept = 1, lin_offset = 1)
  hp_s <- tibble::tibble(
    lin_slope = 1 + 10^(-8),
    lin_intercept = 1,
    lin_offset = 1
  )
  hp_o <- tibble::tibble(
    lin_slope = 1,
    lin_intercept = 1,
    lin_offset = 1 + 10^(-8)
  )

  deriv_s <- lin_kernel(c(1, 2), c(2, 3), hp, "lin_slope") %>% as.vector()
  deriv_o <- lin_kernel(c(1, 2), c(2, 3), hp, "lin_offset")

  emp_deriv_s <- (lin_kernel(c(1, 2), c(2, 3), hp_s)[1] -
    lin_kernel(c(1, 2), c(2, 3), hp)[1]) /
    10^(-8)
  emp_deriv_o <- (lin_kernel(c(1, 2), c(2, 3), hp_o)[1] -
    lin_kernel(c(1, 2), c(2, 3), hp)[1]) /
    10^(-8)

  round(deriv_s, 4) %>% expect_equal(round(emp_deriv_s, 4))
  round(deriv_o, 4) %>% expect_equal(round(emp_deriv_o, 4))
})


## ---- Multi-output convolution kernel ------------------------------------

.mo_inputs <- function() {
  tibble::tibble(
    Output_ID = factor(c(1, 2, 1, 2, 1)),
    Input = c(0.2, 0.5, 0.9, 0.3, 0.7)
  )
}

test_that("convolution_kernel builds a symmetric MO covariance with per-output noise", {
  inp <- .mo_inputs()
  hp <- tibble::tibble(
    Output_ID = factor(1:2), l_t = c(-0.3, 0.2), S_t = c(0.1, -0.2),
    l_u_t = c(0.25, 0.25), noise = c(log(0.5), log(2))
  )
  K <- kern_to_cov(inp, convolution_kernel, hp)
  expect_equal(nrow(K), nrow(inp))
  expect_true(isSymmetric(unname(K), tol = 1e-8))

  K0 <- kern_to_cov(inp, convolution_kernel, dplyr::select(hp, -noise))
  d <- diag(K) - diag(K0)
  expect_equal(unname(d[c(1, 3, 5)]), rep(0.5, 3))  # output 1: exp(log(0.5))
  expect_equal(unname(d[c(2, 4)]), rep(2, 2))        # output 2: exp(log(2))
})

test_that("convolution_kernel single-latent derivatives match finite differences (1D)", {
  inp <- .mo_inputs()
  hp <- tibble::tibble(Output_ID = factor(1:2), l_t = c(-0.3, 0.2),
                       S_t = c(0.1, -0.2), l_u_t = c(0.25, 0.25))
  Kf <- function(h) convolution_kernel(inp, inp, h, vectorized = TRUE)
  fd <- function(col, k) {
    e <- 1e-6; p <- hp; m <- hp
    p[[col]][k] <- p[[col]][k] + e; m[[col]][k] <- m[[col]][k] - e
    (Kf(p) - Kf(m)) / (2 * e)
  }
  for (k in 1:2) {
    expect_equal(convolution_kernel(inp, inp, hp, TRUE, paste0("l_t_", k)),
                 fd("l_t", k), tolerance = 1e-4)
    expect_equal(convolution_kernel(inp, inp, hp, TRUE, paste0("S_t_", k)),
                 fd("S_t", k), tolerance = 1e-4)
  }
  e <- 1e-6; p <- hp; m <- hp
  p$l_u_t <- hp$l_u_t + e; m$l_u_t <- hp$l_u_t - e
  expect_equal(convolution_kernel(inp, inp, hp, TRUE, "l_u_t"),
               (Kf(p) - Kf(m)) / (2 * e), tolerance = 1e-4)
})

test_that("convolution_kernel derivatives match finite differences (2D inputs)", {
  inp <- tibble::tibble(Output_ID = factor(c(1, 2, 1, 2)),
                        Input = c(0.2, 0.9, 0.4, 0.7),
                        Input_2 = c(0.5, 0.1, 0.8, 0.3))
  hp <- tibble::tibble(Output_ID = factor(1:2), l_t = c(-0.3, 0.2),
                       S_t = c(0.1, -0.2), l_u_t = c(0.25, 0.25))
  Kf <- function(h) convolution_kernel(inp, inp, h, vectorized = TRUE)
  fd_lt <- function(k) {
    e <- 1e-6; p <- hp; m <- hp
    p$l_t[k] <- p$l_t[k] + e; m$l_t[k] <- m$l_t[k] - e
    (Kf(p) - Kf(m)) / (2 * e)
  }
  expect_equal(convolution_kernel(inp, inp, hp, TRUE, "l_t_1"), fd_lt(1),
               tolerance = 1e-4)
  e <- 1e-6; p <- hp; m <- hp
  p$l_u_t <- hp$l_u_t + e; m$l_u_t <- hp$l_u_t - e
  expect_equal(convolution_kernel(inp, inp, hp, TRUE, "l_u_t"),
               (Kf(p) - Kf(m)) / (2 * e), tolerance = 1e-4)
})

test_that("convolution_kernel cross-covariance derivatives match finite differences", {
  xa <- tibble::tibble(Output_ID = factor(c(1, 2, 1, 2)), Input = c(0.2, 0.9, 0.4, 0.5))
  xb <- tibble::tibble(Output_ID = factor(c(2, 1, 2)), Input = c(0.3, 0.6, 0.8))
  hp <- tibble::tibble(Output_ID = factor(1:2), l_t = c(-0.3, 0.2),
                       S_t = c(0.1, -0.2), l_u_t = c(0.25, 0.25))
  Kf <- function(h) convolution_kernel(xa, xb, h, vectorized = TRUE)
  fd <- function(col, k) {
    e <- 1e-6; p <- hp; m <- hp
    p[[col]][k] <- p[[col]][k] + e; m[[col]][k] <- m[[col]][k] - e
    (Kf(p) - Kf(m)) / (2 * e)
  }
  ana <- convolution_kernel(xa, xb, hp, TRUE, "S_t_1")
  expect_equal(dim(ana), c(4L, 3L))
  expect_equal(ana, fd("S_t", 1), tolerance = 1e-4)
})

test_that("convolution_kernel multi-latent (Q=2) derivatives match finite differences", {
  inp <- .mo_inputs()
  hp <- tibble::tibble(
    Output_ID = factor(c(1, 2, 1, 2)), Latent_ID = factor(c(1, 1, 2, 2)),
    l_t = c(-0.3, 0.2, 0.1, -0.4), S_t = c(0.1, -0.2, 0.3, 0),
    l_u_t = c(0.25, 0.25, -0.1, -0.1)
  )
  Kf <- function(h) convolution_kernel(inp, inp, h, vectorized = TRUE)
  d_of <- as.integer(as.character(hp$Output_ID))
  q_of <- as.integer(as.character(hp$Latent_ID))
  fd_col <- function(col, d, q) {
    e <- 1e-6; p <- hp; m <- hp; sel <- d_of == d & q_of == q
    p[[col]][sel] <- p[[col]][sel] + e; m[[col]][sel] <- m[[col]][sel] - e
    (Kf(p) - Kf(m)) / (2 * e)
  }
  expect_equal(convolution_kernel(inp, inp, hp, TRUE, "l_t_o1_l2"),
               fd_col("l_t", 1, 2), tolerance = 1e-4)
  expect_equal(convolution_kernel(inp, inp, hp, TRUE, "S_t_o2_l1"),
               fd_col("S_t", 2, 1), tolerance = 1e-4)
  e <- 1e-6; p <- hp; m <- hp; sel <- q_of == 2
  p$l_u_t[sel] <- p$l_u_t[sel] + e; m$l_u_t[sel] <- m$l_u_t[sel] - e
  expect_equal(convolution_kernel(inp, inp, hp, TRUE, "l_u_t_l2"),
               (Kf(p) - Kf(m)) / (2 * e), tolerance = 1e-4)
})

test_that("identify_coord_cols keeps only numeric coordinate columns", {
  df <- tibble::tibble(Task_ID = 1:2, Input_ID = c(1L, 1L), Input = c(0.3, 0.7),
                       Output_ID = factor(c(1, 2)), Output = c(5, 6),
                       Reference = c("a", "b"))
  expect_equal(identify_coord_cols(df), "Input")
})

test_that("flatten/unflatten hp roundtrip (single- and multi-latent)", {
  hp1 <- tibble::tibble(Output_ID = factor(1:2), l_t = c(-0.3, 0.2),
                        S_t = c(0.1, -0.2), l_u_t = c(0.25, 0.25), noise = c(-3, -4))
  f1 <- flatten_hp_mo(hp1)
  expect_setequal(names(f1),
                  c("l_t_1", "l_t_2", "S_t_1", "S_t_2", "l_u_t", "noise_1", "noise_2"))
  b1 <- unflatten_hp_mo(f1)
  expect_equal(as.numeric(b1$l_t), c(-0.3, 0.2))
  expect_equal(unique(b1$l_u_t), 0.25)

  hp2 <- tibble::tibble(
    Output_ID = factor(c(1, 2, 1, 2)), Latent_ID = factor(c(1, 1, 2, 2)),
    l_t = c(-0.3, 0.2, 0.1, -0.4), S_t = c(0.1, -0.2, 0.3, 0),
    l_u_t = c(0.25, 0.25, -0.1, -0.1), noise = c(-3, -4, -3, -4)
  )
  f2 <- flatten_hp_mo(hp2)
  expect_true(all(c("l_t_o1_l1", "l_u_t_l2", "noise_o2") %in% names(f2)))
  b2 <- unflatten_hp_mo(f2)
  expect_true("Latent_ID" %in% names(b2))
  expect_equal(nrow(b2), 4L)
})


test_that("convolution kernel & flatten/unflatten support arbitrary string output labels", {
  inp <- tibble::tibble(Output_ID = c("poids", "taille", "poids"),
                        Input = c(0.2, 0.5, 0.9))
  hp <- tibble::tibble(Output_ID = c("poids", "taille"), l_t = c(-0.3, 0.2),
                       S_t = c(0.1, -0.2), l_u_t = c(0.25, 0.25),
                       noise = c(-3, -4))
  K <- kern_to_cov(inp, convolution_kernel, hp)
  expect_true(isSymmetric(unname(K), tol = 1e-8))

  ## derivative w.r.t. l_t of a string-labelled output matches finite diff
  Kf <- function(h) convolution_kernel(inp, inp, h, vectorized = TRUE)
  e <- 1e-6; p <- hp; m <- hp
  p$l_t[2] <- p$l_t[2] + e; m$l_t[2] <- m$l_t[2] - e
  expect_equal(convolution_kernel(inp, inp, hp, TRUE, "l_t_taille"),
               (Kf(p) - Kf(m)) / (2 * e), tolerance = 1e-4)

  ## flatten uses labels directly; unflatten recovers them as character
  fl <- flatten_hp_mo(hp)
  expect_true(all(c("l_t_poids", "l_t_taille", "l_u_t",
                    "noise_poids", "noise_taille") %in% names(fl)))
  bk <- unflatten_hp_mo(fl)
  expect_type(bk$Output_ID, "character")
  expect_setequal(bk$Output_ID, c("poids", "taille"))
})
