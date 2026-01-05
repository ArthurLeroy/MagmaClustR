# Evidence Lower Bound maximised in MagmaClust

Evidence Lower Bound maximised in MagmaClust

## Usage

``` r
elbo_monitoring_VEM(hp_k, hp_i, db, kern_i, kern_k, hyperpost, m_k, pen_diag)
```

## Arguments

- hp_k:

  A tibble, data frame or named vector of hyper-parameters for each
  clusters.

- hp_i:

  A tibble, data frame or named vector of hyper-parameters for each
  individuals.

- db:

  A tibble containing values we want to compute elbo on. Required
  columns: Input, Output. Additional covariate columns are allowed.

- kern_i:

  Kernel used to compute the covariance matrix of individuals GPs at
  corresponding inputs.

- kern_k:

  Kernel used to compute the covariance matrix of the mean GPs at
  corresponding inputs.

- hyperpost:

  A list of parameters for the variational distributions of the K mean
  GPs.

- m_k:

  Prior value of the mean parameter of the mean GPs (mu_k). Length = 1
  or nrow(db).

- pen_diag:

  A jitter term that is added to the covariance matrix to avoid
  numerical issues when inverting, in cases of nearly singular matrices.

## Value

Value of the elbo that is maximised during the VEM algorithm used for
training in MagmaClust.

## Examples

``` r
TRUE
#> [1] TRUE
```
