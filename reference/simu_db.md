# Simulate a dataset tailored for MagmaClustR

Simulate a complete training dataset, which may be representative of
various applications. Several flexible arguments allow adjustment of the
number of individuals, of observed inputs, and the values of many
parameters controlling the data generation.

## Usage

``` r
simu_db(
  M = 10,
  N = 10,
  K = 1,
  covariate = FALSE,
  grid = seq(0, 10, 0.05),
  grid_cov = seq(0, 10, 0.5),
  common_input = TRUE,
  common_hp = TRUE,
  add_hp = FALSE,
  add_clust = FALSE,
  int_mu_v = c(4, 5),
  int_mu_l = c(0, 1),
  int_i_v = c(1, 2),
  int_i_l = c(0, 1),
  int_i_sigma = c(0, 0.2),
  lambda_int = c(30, 40),
  m_int = c(0, 10),
  lengthscale_int = c(30, 40),
  m0_slope = c(-5, 5),
  m0_intercept = c(-50, 50)
)
```

## Arguments

- M:

  An integer. The number of individual per cluster.

- N:

  An integer. The number of observations per individual.

- K:

  An integer. The number of underlying clusters.

- covariate:

  A logical value indicating whether the dataset should include an
  additional input covariate named 'Covariate'.

- grid:

  A vector of numbers defining a grid of observations (i.e. the
  reference inputs).

- grid_cov:

  A vector of numbers defining a grid of observations (i.e. the
  covariate reference inputs).

- common_input:

  A logical value indicating whether the reference inputs are common to
  all individual.

- common_hp:

  A logical value indicating whether the hyper-parameters are common to
  all individual. If TRUE and K\>1, the hyper-parameters remain
  different between the clusters.

- add_hp:

  A logical value indicating whether the values of hyper-parameters
  should be added as columns in the dataset.

- add_clust:

  A logical value indicating whether the name of the clusters should be
  added as a column in the dataset.

- int_mu_v:

  A vector of 2 numbers, defining an interval of admissible values for
  the variance hyper-parameter of the mean process' kernel.

- int_mu_l:

  A vector of 2 numbers, defining an interval of admissible values for
  the lengthscale hyper-parameter of the mean process' kernel.

- int_i_v:

  A vector of 2 numbers, defining an interval of admissible values for
  the variance hyper-parameter of the individual process' kernel.

- int_i_l:

  A vector of 2 numbers, defining an interval of admissible values for
  the lengthscale hyper-parameter of the individual process' kernel.

- int_i_sigma:

  A vector of 2 numbers, defining an interval of admissible values for
  the noise hyper-parameter.

- lambda_int:

  A vector of 2 numbers, defining an interval of admissible values for
  the lambda parameter of the 2D exponential.

- m_int:

  A vector of 2 numbers, defining an interval of admissible values for
  the mean of the 2D exponential.

- lengthscale_int:

  A vector of 2 numbers, defining an interval of admissible values for
  the lengthscale parameter of the 2D exponential.

- m0_slope:

  A vector of 2 numbers, defining an interval of admissible values for
  the slope of m0.

- m0_intercept:

  A vector of 2 numbers, defining an interval of admissible values for
  the intercept of m0.

## Value

A full dataset of simulated training data.

## Examples

``` r
## Generate a dataset with 3 clusters of 4 individuals, observed at 10 inputs
data = simu_db(M = 4, N = 10, K = 3)

## Generate a 2-D dataset with an additional input 'Covariate'
data = simu_db(covariate = TRUE)

## Generate a dataset where input locations are different among individuals
data = simu_db(common_input = FALSE)

## Generate a dataset with an additional column indicating the true clusters
data = simu_db(K = 3, add_clust = TRUE)
```
