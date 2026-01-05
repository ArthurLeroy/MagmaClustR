# Modified log-Likelihood function for GPs

Log-Likelihood function involved in Magma during the maximisation step
of the training. The log-Likelihood is defined as a simple Gaussian
likelihood added with correction trace term.

## Usage

``` r
logL_GP_mod(hp, db, mean, kern, post_cov, pen_diag)
```

## Arguments

- hp:

  A tibble, data frame or named vector of hyper-parameters.

- db:

  A tibble containing values we want to compute logL on. Required
  columns: Input, Output. Additional covariate columns are allowed.

- mean:

  A vector, specifying the mean of the GP at the reference inputs.

- kern:

  A kernel function.

- post_cov:

  A matrix, covariance parameter of the hyper-posterior. Used to compute
  the correction term.

- pen_diag:

  A jitter term that is added to the covariance matrix to avoid
  numerical issues when inverting, in cases of nearly singular matrices.

## Value

A number, corresponding to the value of the modified Gaussian
log-Likelihood defined in Magma.

## Examples

``` r
TRUE
#> [1] TRUE
```
