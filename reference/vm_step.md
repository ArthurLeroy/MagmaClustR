# V-Step of the VEM algorithm

Maximization step of the Variational EM algorithm used to compute
hyper-parameters of all the kernels involved in MagmaClust.

## Usage

``` r
vm_step(
  db,
  old_hp_k,
  old_hp_i,
  list_mu_param,
  kern_k,
  kern_i,
  m_k,
  common_hp_k,
  common_hp_i,
  pen_diag
)
```

## Arguments

- db:

  A tibble or data frame. Columns required: ID, Input, Output.
  Additional columns for covariates can be specified.

- old_hp_k:

  A named vector, tibble or data frame, containing the hyper-parameters
  from the previous M-step (or initialisation) associated with the mean
  GPs.

- old_hp_i:

  A named vector, tibble or data frame, containing the hyper-parameters
  from the previous M-step (or initialisation) associated with the
  individual GPs.

- list_mu_param:

  List of parameters of the K mean GPs.

- kern_k:

  A kernel used to compute the covariance matrix of the mean GP at
  corresponding timestamps.

- kern_i:

  A kernel used to compute the covariance matrix of individuals GP at
  corresponding timestamps.

- m_k:

  A named list of prior mean parameters for the K mean GPs. Length = 1
  or nrow(unique(db\$Input))

- common_hp_k:

  A boolean indicating whether hp are common among mean GPs (for each
  mu_k)

- common_hp_i:

  A boolean indicating whether hp are common among individual GPs (for
  each y_i)

- pen_diag:

  A number. A jitter term, added on the diagonal to prevent numerical
  issues when inverting nearly singular matrices.

## Value

A named list, containing the elements `hp_k`, a tibble containing the
hyper-parameters associated with each cluster, `hp_i`, a tibble
containing the hyper-parameters associated with the individual GPs, and
`prop_mixture_k`, a tibble containing the hyper-parameters associated
with each individual, indicating the probabilities to belong to each
cluster.

## Examples

``` r
TRUE
#> [1] TRUE
```
