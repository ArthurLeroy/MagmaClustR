# Penalised elbo for multiple mean GPs with common HPs

Penalised elbo for multiple mean GPs with common HPs

## Usage

``` r
elbo_GP_mod_common_hp_k(hp, db, mean, kern, post_cov, pen_diag)
```

## Arguments

- hp:

  A tibble, data frame or named vector containing hyper-parameters.

- db:

  A tibble containing values we want to compute elbo on. Required
  columns: Input, Output. Additional covariate columns are allowed.

- mean:

  A list of the K mean GPs at union of observed timestamps.

- kern:

  A kernel function used to compute the covariance matrix at
  corresponding timestamps.

- post_cov:

  A List of the K posterior covariance of the mean GP (mu_k). Used to
  compute correction term (cor_term).

- pen_diag:

  A jitter term that is added to the covariance matrix to avoid
  numerical issues when inverting, in cases of nearly singular matrices.

## Value

The value of the penalised Gaussian elbo for the sum of the k mean GPs
with common HPs.

## Examples

``` r
TRUE
#> [1] TRUE
```
