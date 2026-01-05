# E-Step of the EM algorithm

Expectation step of the EM algorithm to compute the parameters of the
hyper-posterior Gaussian distribution of the mean process in Magma.

## Usage

``` r
e_step(db, m_0, kern_0, kern_i, hp_0, hp_i, pen_diag)
```

## Arguments

- db:

  A tibble or data frame. Columns required: ID, Input, Output.
  Additional columns for covariates can be specified.

- m_0:

  A vector, corresponding to the prior mean of the mean GP.

- kern_0:

  A kernel function, associated with the mean GP.

- kern_i:

  A kernel function, associated with the individual GPs.

- hp_0:

  A named vector, tibble or data frame of hyper-parameters associated
  with `kern_0`.

- hp_i:

  A tibble or data frame of hyper-parameters associated with `kern_i`.

- pen_diag:

  A number. A jitter term, added on the diagonal to prevent numerical
  issues when inverting nearly singular matrices.

## Value

A named list, containing the elements `mean`, a tibble containing the
Input and associated Output of the hyper-posterior's mean parameter, and
`cov`, the hyper-posterior's covariance matrix.

## Examples

``` r
TRUE
#> [1] TRUE
```
