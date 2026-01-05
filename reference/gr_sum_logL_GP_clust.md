# Gradient of the mixture of Gaussian likelihoods

Compute the gradient of a sum of Gaussian log-likelihoods, weighted by
their mixture probabilities.

## Usage

``` r
gr_sum_logL_GP_clust(hp, db, mixture, mean, kern, post_cov, pen_diag)
```

## Arguments

- hp:

  A tibble, data frame or named vector of hyper-parameters.

- db:

  A tibble containing data we want to evaluate the logL on. Required
  columns: Input, Output. Additional covariate columns are allowed.

- mixture:

  A tibble or data frame, indicating the mixture probabilities of each
  cluster for the new individual/task.

- mean:

  A list of hyper-posterior mean parameters for all clusters.

- kern:

  A kernel function.

- post_cov:

  A list of hyper-posterior covariance parameters for all clusters.

- pen_diag:

  A jitter term that is added to the covariance matrix to avoid
  numerical issues when inverting, in cases of nearly singular matrices.

## Value

A named vector, corresponding to the value of the hyper-parameters'
gradients for the mixture of Gaussian log-likelihoods involved in the
prediction step of MagmaClust.

## Examples

``` r
TRUE
#> [1] TRUE
```
