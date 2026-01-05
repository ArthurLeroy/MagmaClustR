# Gradient of the modified logLikelihood with common HPs for GPs in Magma

Gradient of the modified logLikelihood with common HPs for GPs in Magma

## Usage

``` r
gr_GP_mod_common_hp(hp, db, mean, kern, post_cov, pen_diag)
```

## Arguments

- hp:

  A tibble or data frame containing hyper-parameters for all
  individuals.

- db:

  A tibble containing the values we want to compute the logL on.
  Required columns: ID, Input, Output. Additional covariate columns are
  allowed.

- mean:

  A vector, specifying the mean of the GPs at the reference inputs.

- kern:

  A kernel function.

- post_cov:

  A matrix, covariance parameter of the hyper-posterior. Used to compute
  the correction term.

- pen_diag:

  A jitter term that is added to the covariance matrix to avoid
  numerical issues when inverting, in cases of nearly singular matrices.

## Value

A named vector, corresponding to the value of the hyper-parameters'
gradients for the modified Gaussian log-Likelihood involved in Magma
with the 'common HP' setting.

## Examples

``` r
TRUE
#> [1] TRUE
```
