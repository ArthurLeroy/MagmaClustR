# Gradient of the elbo for a mixture of GPs

Gradient of the elbo for a mixture of GPs

## Usage

``` r
gr_clust_multi_GP(hp, db, hyperpost, kern, pen_diag)
```

## Arguments

- hp:

  A tibble, data frame or named vector containing hyper-parameters.

- db:

  A tibble containing the values we want to compute the elbo on.
  Required columns: Input, Output. Additional covariate columns are
  allowed.

- hyperpost:

  List of parameters for the K mean Gaussian processes.

- kern:

  A kernel function.

- pen_diag:

  A jitter term that is added to the covariance matrix to avoid
  numerical issues when inverting, in cases of nearly singular matrices.

## Value

The gradient of the penalised Gaussian elbo for a mixture of GPs

## Examples

``` r
TRUE
#> [1] TRUE
```
