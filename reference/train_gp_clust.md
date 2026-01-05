# Prediction in MagmaClust: learning new HPs and mixture probabilities

Learning hyper-parameters and mixture probabilities of any new
individual/task is required in `MagmaClust` in the prediction procedure.
By providing data for the new individual/task, the hyper-posterior mean
and covariance parameters, the mixture proportions, and initialisation
values for the hyper-parameters, `train_gp_clust` uses an EM algorithm
to compute maximum likelihood estimates of the hyper-parameters and
hyper-posterior mixture probabilities of the new individual/task.

## Usage

``` r
train_gp_clust(
  data,
  prop_mixture = NULL,
  ini_hp = NULL,
  kern = "SE",
  hyperpost = NULL,
  pen_diag = 1e-10,
  n_iter_max = 25,
  cv_threshold = 0.001
)
```

## Arguments

- data:

  A tibble or data frame. Required columns: `Input`, `Output`.
  Additional columns for covariates can be specified. The `Input` column
  should define the variable that is used as reference for the
  observations (e.g. time for longitudinal data). The `Output` column
  specifies the observed values (the response variable). The data frame
  can also provide as many covariates as desired, with no constraints on
  the column names. These covariates are additional inputs (explanatory
  variables) of the models that are also observed at each reference
  `Input`.

- prop_mixture:

  A tibble or a named vector. Each name of column or element should
  refer to a cluster. The value associated with each cluster is a number
  between 0 and 1, corresponding to the mixture proportions.

- ini_hp:

  A tibble or data frame of hyper-parameters associated with `kern`, the
  individual process kernel. The
  [`hp`](https://arthurleroy.github.io/MagmaClustR/reference/hp.md)
  function can be used to draw custom hyper-parameters with the correct
  format.

- kern:

  A kernel function, defining the covariance structure of the GP.
  Several popular kernels (see [The Kernel
  Cookbook](https://www.cs.toronto.edu/~duvenaud/cookbook/)) are already
  implemented and can be selected within the following list:

  - "SE": (default value) the Squared Exponential Kernel (also called
    Radial Basis Function or Gaussian kernel),

  - "LIN": the Linear kernel,

  - "PERIO": the Periodic kernel,

  - "RQ": the Rational Quadratic kernel. Compound kernels can be created
    as sums or products of the above kernels. For combining kernels,
    simply provide a formula as a character string where elements are
    separated by whitespaces (e.g. "SE + PERIO"). As the elements are
    treated sequentially from the left to the right, the product
    operator '\*' shall always be used before the '+' operators (e.g.
    'SE \* LIN + RQ' is valid whereas 'RQ + SE \* LIN' is not).

- hyperpost:

  A list, containing the elements `mean`, `cov` and `mixture` the
  parameters of the hyper-posterior distributions of the mean processes.
  Typically, this argument should come from a previous learning using
  [`train_magmaclust`](https://arthurleroy.github.io/MagmaClustR/reference/train_magmaclust.md),
  or a previous prediction with
  [`pred_magmaclust`](https://arthurleroy.github.io/MagmaClustR/reference/pred_magmaclust.md),
  with the argument `get_hyperpost` set to TRUE.

- pen_diag:

  A number. A jitter term, added on the diagonal to prevent numerical
  issues when inverting nearly singular matrices.

- n_iter_max:

  A number, indicating the maximum number of iterations of the EM
  algorithm to proceed while not reaching convergence.

- cv_threshold:

  A number, indicating the threshold of the likelihood gain under which
  the EM algorithm will stop.

## Value

A list, containing the results of the EM algorithm used during the
prediction step of MagmaClust. The elements of the list are:

- hp: A tibble of optimal hyper-parameters for the new individual's GP.

- mixture: A tibble of mixture probabilities for the new individual.

## Examples

``` r
TRUE
#> [1] TRUE
```
