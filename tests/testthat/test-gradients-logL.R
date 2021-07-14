test_that("gradient of logL_GP_mod() works for Squared Exponential kernel", {
  db = tibble::tibble(Input= 1:5, Output= 2:6)
  mean = rep(0, 5)
  hp = tibble::tibble(variance = 1, lengthscale = 0.5)
  new_cov = kern_to_cov(db$Input, 'SE', hp)

  hp_v = tibble::tibble(variance = 1 + 10^(-8), lengthscale = 0.5)
  hp_l = tibble::tibble(variance = 1, lengthscale = 0.5 + 10^(-8))

  deriv_v = gr_GP_mod(hp, db, mean, 'SE', new_cov, 0)[['variance']]
  deriv_l = gr_GP_mod(hp, db, mean, 'SE', new_cov, 0)[['lengthscale']]

  emp_deriv_v = (logL_GP_mod(hp_v, db, mean, 'SE', new_cov, 0) -
               logL_GP_mod(hp, db, mean, 'SE', new_cov, 0)) / 10^(-8)
  emp_deriv_l = (logL_GP_mod(hp_l, db, mean, 'SE', new_cov, 0) -
                logL_GP_mod(hp, db, mean, 'SE', new_cov, 0)) / 10^(-8)

  round(deriv_v, 4) %>% expect_equal(round(emp_deriv_v, 4))
  round(deriv_l, 4) %>% expect_equal(round(emp_deriv_l, 4))
})
