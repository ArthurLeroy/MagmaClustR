# Modified log-Likelihood function with common HPs for GPs

Log-Likelihood function involved in Magma during the maximisation step
of the training, in the particular case where the hyper-parameters are
shared by all individuals. The log-Likelihood is defined as a sum over
all individuals of Gaussian likelihoods added with correction trace
terms.

## Usage

``` r
logL_GP_mod_common_hp(hp, db, mean, kern, post_cov, pen_diag)
```

## Arguments

- hp:

  A tibble, data frame of hyper-parameters.

- db:

  A tibble containing the values we want to compute the logL on.
  Required columns: ID, Input, Output. Additional covariate columns are
  allowed.

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
log-Likelihood with common hyper-parameters defined in Magma.

## Examples

``` r
TRUE
#> [1] TRUE
```
