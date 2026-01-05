# Compute a mixture of Gaussian log-likelihoods

During the prediction step of MagmaClust, an EM algorithm is used to
compute the maximum likelihood estimator of the hyper-parameters along
with mixture probabilities for the new individual/task. This function
implements the quantity that is maximised (i.e. a sum of Gaussian
log-likelihoods, weighted by their mixture probabilities). It can also
be used to monitor the EM algorithm when providing the 'prop_mixture'
argument, for proper penalisation of the full log-likelihood.

## Usage

``` r
sum_logL_GP_clust(
  hp,
  db,
  mixture,
  mean,
  kern,
  post_cov,
  prop_mixture = NULL,
  pen_diag
)
```

## Arguments

- hp:

  A tibble, data frame or named vector of hyper-parameters.

- db:

  A tibble containing data we want to evaluate the logL on. Required
  columns: Input, Output. Additional covariate columns are allowed.

- mixture:

  A tibble or data frame, indicating the mixture probabilities of each
  cluster for the new individual/task.

- mean:

  A list of hyper-posterior mean parameters for all clusters.

- kern:

  A kernel function.

- post_cov:

  A list of hyper-posterior covariance parameters for all clusters.

- prop_mixture:

  A tibble or a named vector. Each name of column or element should
  refer to a cluster. The value associated with each cluster is a number
  between 0 and 1, corresponding to the mixture proportions.

- pen_diag:

  A jitter term that is added to the covariance matrix to avoid
  numerical issues when inverting, in cases of nearly singular matrices.

## Value

A number, expectation of mixture of Gaussian log-likelihoods in the
prediction step of MagmaClust. This quantity is supposed to increase at
each step of the EM algorithm, and can be used for monitoring the
procedure.

## Examples

``` r
TRUE
#> [1] TRUE
```
