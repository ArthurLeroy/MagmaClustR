# Update the mixture probabilities for each individual and each cluster

Update the mixture probabilities for each individual and each cluster

## Usage

``` r
update_mixture(db, mean_k, cov_k, hp, kern, prop_mixture, pen_diag)
```

## Arguments

- db:

  A tibble or data frame. Columns required: `ID`, `Input`, `Output`.
  Additional columns for covariates can be specified.

- mean_k:

  A list of the K hyper-posterior mean parameters.

- cov_k:

  A list of the K hyper-posterior covariance matrices.

- hp:

  A named vector, tibble or data frame of hyper-parameters associated
  with `kern`, the individual process' kernel. The columns/elements
  should be named according to the hyper-parameters that are used in
  `kern`.

- kern:

  A kernel function, defining the covariance structure of the individual
  GPs.

- prop_mixture:

  A tibble containing the hyper-parameters associated with each
  individual, indicating in which cluster it belongs.

- pen_diag:

  A number. A jitter term, added on the diagonal to prevent numerical
  issues when inverting nearly singular matrices.

## Value

Compute the hyper-posterior multinomial distributions by updating
mixture probabilities.

## Examples

``` r
TRUE
#> [1] TRUE
```
