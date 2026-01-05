# Penalised elbo for multiple individual GPs with common HPs

Penalised elbo for multiple individual GPs with common HPs

## Usage

``` r
elbo_clust_multi_GP_common_hp_i(hp, db, hyperpost, kern, pen_diag)
```

## Arguments

- hp:

  A tibble, data frame or named vector containing hyper-parameters.

- db:

  A tibble containing values we want to compute elbo on. Required
  columns: Input, Output. Additional covariate columns are allowed.

- hyperpost:

  List of parameters for the K mean Gaussian processes.

- kern:

  A kernel function used to compute the covariance matrix at
  corresponding timestamps.

- pen_diag:

  A jitter term that is added to the covariance matrix to avoid
  numerical issues when inverting, in cases of nearly singular matrices.

## Value

The value of the penalised Gaussian elbo for the sum of the M individual
GPs with common HPs.

## Examples

``` r
TRUE
#> [1] TRUE
```
