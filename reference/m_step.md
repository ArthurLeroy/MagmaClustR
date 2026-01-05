# M-Step of the EM algorithm

Maximisation step of the EM algorithm to compute hyper-parameters of all
the kernels involved in Magma.

## Usage

``` r
m_step(
  db,
  m_0,
  kern_0,
  kern_i,
  old_hp_0,
  old_hp_i,
  post_mean,
  post_cov,
  common_hp,
  pen_diag
)
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

- old_hp_0:

  A named vector, tibble or data frame, containing the hyper-parameters
  from the previous M-step (or initialisation) associated with the mean
  GP.

- old_hp_i:

  A tibble or data frame, containing the hyper-parameters from the
  previous M-step (or initialisation) associated with the individual
  GPs.

- post_mean:

  A tibble, coming out of the E step, containing the Input and
  associated Output of the hyper-posterior mean parameter.

- post_cov:

  A matrix, coming out of the E step, being the hyper-posterior
  covariance parameter.

- common_hp:

  A logical value, indicating whether the set of hyper-parameters is
  assumed to be common to all indiviuals.

- pen_diag:

  A number. A jitter term, added on the diagonal to prevent numerical
  issues when inverting nearly singular matrices.

## Value

A named list, containing the elements `hp_0`, a tibble containing the
hyper-parameters associated with the mean GP, `hp_i`, a tibble
containing the hyper-parameters associated with the individual GPs.

## Examples

``` r
TRUE
#> [1] TRUE
```
