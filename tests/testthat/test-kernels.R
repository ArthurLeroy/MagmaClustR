test_that("Squared Exponential kernel works for null distance", {
  hp <- tibble::tibble(variance = 1, lengthscale = 1)

  se_kernel(2, 2, hp)[1] %>% expect_equal(exp(1))
})

test_that("Periodic kernel works for null distance", {
  hp <- tibble::tibble(variance = 1, lengthscale = 1, period = 1)

  perio_kernel(2, 2, hp)[1] %>% expect_equal(exp(1))
})

test_that("Rational quadratic kernel works for null distance", {
  hp <- tibble::tibble(variance = 1, lengthscale = 1, scale = 1)

  rq_kernel(2, 2, hp)[1] %>% expect_equal(exp(1))
})


test_that("gradients for the Squared Exponential kernel are valid", {
  hp <- tibble::tibble(variance = 1, lengthscale = 1)
  hp_v <- tibble::tibble(variance = 1 + 10^(-8), lengthscale = 1)
  hp_l <- tibble::tibble(variance = 1, lengthscale = 1 + 10^(-8))

  deriv_v = se_kernel(c(1,2), c(2,3), hp, 'variance')
  deriv_l = se_kernel(c(1,2), c(2,3), hp, 'lengthscale')

  emp_deriv_v = (se_kernel(c(1,2), c(2,3), hp_v)[1] -
                   se_kernel(c(1,2), c(2,3), hp)[1]) / 10^(-8)
  emp_deriv_l = (se_kernel(c(1,2), c(2,3), hp_l)[1] -
                   se_kernel(c(1,2), c(2,3), hp)[1]) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(emp_deriv_v, 4))
  round(deriv_l, 4) %>% expect_equal(round(emp_deriv_l, 4))
})

test_that("gradients for the Periodic kernel are valid", {
  hp <- tibble::tibble(variance = 1, lengthscale = 1, period = pi)
  hp_v <- tibble::tibble(variance = 1 + 10^(-8), lengthscale = 1, period = pi)
  hp_l <- tibble::tibble(variance = 1, lengthscale = 1 + 10^(-8), period = pi)
  hp_p <- tibble::tibble(variance = 1, lengthscale = 1, period = pi + 10^(-8))

  deriv_v = perio_kernel(c(1,2), c(2,3), hp, 'variance')
  deriv_l = perio_kernel(c(1,2), c(2,3), hp, 'lengthscale')
  deriv_p = perio_kernel(c(1,2), c(2,3), hp, 'period')

  emp_deriv_v = (perio_kernel(c(1,2), c(2,3), hp_v)[1] -
                   perio_kernel(c(1,2), c(2,3), hp)[1]) / 10^(-8)
  emp_deriv_l = (perio_kernel(c(1,2), c(2,3), hp_l)[1] -
                   perio_kernel(c(1,2), c(2,3), hp)[1]) / 10^(-8)
  emp_deriv_p = (perio_kernel(c(1,2), c(2,3), hp_p)[1] -
                   perio_kernel(c(1,2), c(2,3), hp)[1]) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(emp_deriv_v, 4))
  round(deriv_l, 4) %>% expect_equal(round(emp_deriv_l, 4))
  round(deriv_p, 4) %>% expect_equal(round(emp_deriv_p, 4))
})

test_that("gradients for the Rational Quadratic kernel are valid", {
  hp <- tibble::tibble(variance = 1, lengthscale = 1, scale = 1)
  hp_v <- tibble::tibble(variance = 1 + 10^(-8), lengthscale = 1, scale = 1)
  hp_l <- tibble::tibble(variance = 1, lengthscale = 1 + 10^(-8), scale = 1)
  hp_s <- tibble::tibble(variance = 1, lengthscale = 1, scale = 1 + 10^(-8))

  deriv_v = rq_kernel(c(1,2), c(2,3), hp, 'variance')
  deriv_l = rq_kernel(c(1,2), c(2,3), hp, 'lengthscale')
  deriv_s = rq_kernel(c(1,2), c(2,3), hp, 'scale')

  emp_deriv_v = (rq_kernel(c(1,2), c(2,3), hp_v)[1] -
                   rq_kernel(c(1,2), c(2,3), hp)[1]) / 10^(-8)
  emp_deriv_l = (rq_kernel(c(1,2), c(2,3), hp_l)[1] -
                   rq_kernel(c(1,2), c(2,3), hp)[1]) / 10^(-8)
  emp_deriv_s = (rq_kernel(c(1,2), c(2,3), hp_s)[1] -
                   rq_kernel(c(1,2), c(2,3), hp)[1]) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(emp_deriv_v, 4))
  round(deriv_l, 4) %>% expect_equal(round(emp_deriv_l, 4))
  round(deriv_s, 4) %>% expect_equal(round(emp_deriv_s, 4))
})
