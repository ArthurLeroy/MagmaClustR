# Log-Likelihood for monitoring the EM algorithm in Magma

Log-Likelihood for monitoring the EM algorithm in Magma

## Usage

``` r
logL_monitoring(
  hp_0,
  hp_i,
  db,
  m_0,
  kern_0,
  kern_i,
  post_mean,
  post_cov,
  pen_diag
)
```

## Arguments

- hp_0:

  A named vector, tibble or data frame, containing the hyper-parameters
  associated with the mean GP.

- hp_i:

  A tibble or data frame, containing the hyper-parameters with the
  individual GPs.

- db:

  A tibble or data frame. Columns required: ID, Input, Output.
  Additional columns for covariates can be specified.

- m_0:

  A vector, corresponding to the prior mean of the mean GP.

- kern_0:

  A kernel function, associated with the mean GP.

- kern_i:

  A kernel function, associated with the individual GPs.

- post_mean:

  A tibble, coming out of the E step, containing the Input and
  associated Output of the hyper-posterior mean parameter.

- post_cov:

  A matrix, coming out of the E step, being the hyper-posterior
  covariance parameter.

- pen_diag:

  A jitter term that is added to the covariance matrix to avoid
  numerical issues when inverting, in cases of nearly singular matrices.

## Value

A number, expectation of joint log-likelihood of the model. This
quantity is supposed to increase at each step of the EM algorithm, and
thus used for monitoring the procedure.

## Examples

``` r
TRUE
#> [1] TRUE
```
