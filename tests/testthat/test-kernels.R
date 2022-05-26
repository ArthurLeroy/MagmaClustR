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
    se_kernel(c(1, 2), c(2, 3), hp)[1]) / 10^(-8)
  emp_deriv_l <- (se_kernel(c(1, 2), c(2, 3), hp_l)[1] -
    se_kernel(c(1, 2), c(2, 3), hp)[1]) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(emp_deriv_v, 4))
  round(deriv_l, 4) %>% expect_equal(round(emp_deriv_l, 4))
})

test_that("gradients for the Periodic kernel are valid", {
  hp <- tibble::tibble(
    perio_variance = 1,
    perio_lengthscale = 1, period = pi
  )
  hp_v <- tibble::tibble(
    perio_variance = 1 + 10^(-8),
    perio_lengthscale = 1, period = pi
  )
  hp_l <- tibble::tibble(
    perio_variance = 1,
    perio_lengthscale = 1 + 10^(-8), period = pi
  )
  hp_p <- tibble::tibble(
    perio_variance = 1,
    perio_lengthscale = 1, period = pi + 10^(-8)
  )

  deriv_v <- perio_kernel(c(1, 2), c(2, 3), hp, "perio_variance")
  deriv_l <- perio_kernel(c(1, 2), c(2, 3), hp, "perio_lengthscale")
  deriv_p <- perio_kernel(c(1, 2), c(2, 3), hp, "period")

  emp_deriv_v <- (perio_kernel(c(1, 2), c(2, 3), hp_v)[1] -
    perio_kernel(c(1, 2), c(2, 3), hp)[1]) / 10^(-8)
  emp_deriv_l <- (perio_kernel(c(1, 2), c(2, 3), hp_l)[1] -
    perio_kernel(c(1, 2), c(2, 3), hp)[1]) / 10^(-8)
  emp_deriv_p <- (perio_kernel(c(1, 2), c(2, 3), hp_p)[1] -
    perio_kernel(c(1, 2), c(2, 3), hp)[1]) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(emp_deriv_v, 4))
  round(deriv_l, 4) %>% expect_equal(round(emp_deriv_l, 4))
  round(deriv_p, 4) %>% expect_equal(round(emp_deriv_p, 4))
})

test_that("gradients for the Rational Quadratic kernel are valid", {
  hp <- tibble::tibble(rq_variance = 1, rq_lengthscale = 1, rq_scale = 1)
  hp_v <- tibble::tibble(
    rq_variance = 1 + 10^(-8),
    rq_lengthscale = 1, rq_scale = 1
  )
  hp_l <- tibble::tibble(
    rq_variance = 1,
    rq_lengthscale = 1 + 10^(-8), rq_scale = 1
  )
  hp_s <- tibble::tibble(
    rq_variance = 1,
    rq_lengthscale = 1, rq_scale = 1 + 10^(-8)
  )

  deriv_v <- rq_kernel(c(1, 2), c(2, 3), hp, "rq_variance")
  deriv_l <- rq_kernel(c(1, 2), c(2, 3), hp, "rq_lengthscale")
  deriv_s <- rq_kernel(c(1, 2), c(2, 3), hp, "rq_scale")

  emp_deriv_v <- (rq_kernel(c(1, 2), c(2, 3), hp_v)[1] -
    rq_kernel(c(1, 2), c(2, 3), hp)[1]) / 10^(-8)
  emp_deriv_l <- (rq_kernel(c(1, 2), c(2, 3), hp_l)[1] -
    rq_kernel(c(1, 2), c(2, 3), hp)[1]) / 10^(-8)
  emp_deriv_s <- (rq_kernel(c(1, 2), c(2, 3), hp_s)[1] -
    rq_kernel(c(1, 2), c(2, 3), hp)[1]) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(emp_deriv_v, 4))
  round(deriv_l, 4) %>% expect_equal(round(emp_deriv_l, 4))
  round(deriv_s, 4) %>% expect_equal(round(emp_deriv_s, 4))
})

test_that("gradients for the Linear kernel are valid", {
  hp <- tibble::tibble(lin_slope = 1, lin_intercept = 1, lin_offset = 1)
  hp_s <- tibble::tibble(
    lin_slope = 1 + 10^(-8),
    lin_intercept = 1, lin_offset = 1
  )
  hp_o <- tibble::tibble(
    lin_slope = 1,
    lin_intercept = 1, lin_offset = 1 + 10^(-8)
  )

  deriv_s <- lin_kernel(c(1, 2), c(2, 3), hp, "lin_slope") %>% as.vector()
  deriv_o <- lin_kernel(c(1, 2), c(2, 3), hp, "lin_offset")

  emp_deriv_s <- (lin_kernel(c(1, 2), c(2, 3), hp_s)[1] -
    lin_kernel(c(1, 2), c(2, 3), hp)[1]) / 10^(-8)
  emp_deriv_o <- (lin_kernel(c(1, 2), c(2, 3), hp_o)[1] -
    lin_kernel(c(1, 2), c(2, 3), hp)[1]) / 10^(-8)

  round(deriv_s, 4) %>% expect_equal(round(emp_deriv_s, 4))
  round(deriv_o, 4) %>% expect_equal(round(emp_deriv_o, 4))
})
