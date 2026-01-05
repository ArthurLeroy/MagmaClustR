# Log-Likelihood function of a Gaussian Process

Log-Likelihood function of a Gaussian Process

## Usage

``` r
logL_GP(hp, db, mean, kern, post_cov, pen_diag)
```

## Arguments

- hp:

  A tibble, data frame or named vector containing hyper-parameters.

- db:

  A tibble containing the values we want to compute the logL on.
  Required columns: Input, Output. Additional covariate columns are
  allowed.

- mean:

  A vector, specifying the mean of the GP at the reference inputs.

- kern:

  A kernel function.

- post_cov:

  (optional) A matrix, corresponding to covariance parameter of the
  hyper-posterior. Used to compute the hyper-prior distribution of a new
  individual in Magma.

- pen_diag:

  A jitter term that is added to the covariance matrix to avoid
  numerical issues when inverting, in cases of nearly singular matrices.

## Value

A number, corresponding to the value of Gaussian log-Likelihood (where
the covariance can be the sum of the individual and the
hyper-posterior's mean process covariances).

## Examples

``` r
TRUE
#> [1] TRUE
```
