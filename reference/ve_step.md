# E-Step of the VEM algorithm

Expectation step of the Variational EM algorithm used to compute the
parameters of the hyper-posteriors distributions for the mean processes
and mixture variables involved in MagmaClust.

## Usage

``` r
ve_step(db, m_k, kern_k, kern_i, hp_k, hp_i, old_mixture, iter, pen_diag)
```

## Arguments

- db:

  A tibble or data frame. Columns required: ID, Input, Output.
  Additional columns for covariates can be specified.

- m_k:

  A named list of vectors, corresponding to the prior mean parameters of
  the K mean GPs.

- kern_k:

  A kernel function, associated with the K mean GPs.

- kern_i:

  A kernel function, associated with the M individual GPs.

- hp_k:

  A named vector, tibble or data frame of hyper-parameters associated
  with `kern_k`.

- hp_i:

  A named vector, tibble or data frame of hyper-parameters associated
  with `kern_i`.

- old_mixture:

  A list of mixture values from the previous iteration.

- iter:

  A number, indicating the current iteration of the VEM algorithm.

- pen_diag:

  A number. A jitter term, added on the diagonal to prevent numerical
  issues when inverting nearly singular matrices.

## Value

A named list, containing the elements `mean`, a tibble containing the
Input and associated Output of the hyper-posterior mean parameters,
`cov`, the hyper-posterior covariance matrices, and `mixture`, the
probabilities to belong to each cluster for each individual.

## Examples

``` r
TRUE
#> [1] TRUE
```
