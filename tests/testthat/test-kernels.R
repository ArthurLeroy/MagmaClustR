test_that("Squared Exponential kernel works for null distance", {
  hp <- tibble::tibble(variance = 1, lengthscale = 1)
  expect_equal(se_kernel(2, 2, hp)[1], exp(1))
})

test_that("Periodic kernel works for null distance", {
  hp <- tibble::tibble(variance = 1, lengthscale = 1, period = 1)
  expect_equal(perio_kernel(2, 2, hp)[1], exp(1))
})

test_that("Rational quadratic kernel works for null distance", {
  hp <- tibble::tibble(variance = 1, lengthscale = 1, scale = 1)
  expect_equal(rq_kernel(2, 2, hp)[1], exp(1))
})
