# Compute the hyper-posterior distribution for each cluster in MagmaClust

Recompute the E-step of the VEM algorithm in MagmaClust for a new set of
reference `Input`. Once training is completed, it can be necessary to
evaluate the hyper-posterior distributions of the mean processes at
specific locations, for which we want to make predictions. This process
is directly implemented in the
[`pred_magmaclust`](https://arthurleroy.github.io/MagmaClustR/reference/pred_magmaclust.md)
function but the user might want to use `hyperpost_clust` for a tailored
control of the prediction procedure.

## Usage

``` r
hyperposterior_clust(
  trained_model = NULL,
  data = NULL,
  mixture = NULL,
  hp_k = NULL,
  hp_i = NULL,
  kern_k = NULL,
  kern_i = NULL,
  prior_mean_k = NULL,
  grid_inputs = NULL,
  pen_diag = 1e-10
)
```

## Arguments

- trained_model:

  A list, containing the information coming from a Magma model,
  previously trained using the
  [`train_magma`](https://arthurleroy.github.io/MagmaClustR/reference/train_magma.md)
  function. If `trained_model` is not provided, the arguments `data`,
  `mixture`, `hp_k`, `hp_i`, `kern_k`, and `kern_i` are all required.

- data:

  A tibble or data frame. Required columns: `ID`, `Input` , `Output`.
  Additional columns for covariates can be specified. The `ID` column
  contains the unique names/codes used to identify each individual/task
  (or batch of data). The `Input` column should define the variable that
  is used as reference for the observations (e.g. time for longitudinal
  data). The `Output` column specifies the observed values (the response
  variable). The data frame can also provide as many covariates as
  desired, with no constraints on the column names. These covariates are
  additional inputs (explanatory variables) of the models that are also
  observed at each reference `Input`. Recovered from `trained_model` if
  not provided.

- mixture:

  A tibble or data frame, indicating the mixture probabilities of each
  cluster for each individual. Required column: `ID`. Recovered from
  `trained_model` if not provided.

- hp_k:

  A tibble or data frame of hyper-parameters associated with `kern_k`.
  Recovered from `trained_model` if not provided.

- hp_i:

  A tibble or data frame of hyper-parameters associated with `kern_i`.
  Recovered from `trained_model` if not provided.

- kern_k:

  A kernel function, associated with the mean GPs. Several popular
  kernels (see [The Kernel
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
    Recovered from `trained_model` if not provided.

- kern_i:

  A kernel function, associated with the individual GPs. ("SE", "LIN",
  PERIO" and "RQ" are also available here). Recovered from
  `trained_model` if not provided.

- prior_mean_k:

  The set of hyper-prior mean parameters (m_k) for the K mean GPs, one
  value for each cluster. cluster. This argument can be specified under
  various formats, such as:

  - NULL (default). All hyper-prior means would be set to 0 everywhere.

  - A numerical vector of the same length as the number of clusters.
    Each number is associated with one cluster, and considered to be the
    hyper-prior mean parameter of the cluster (i.e. a constant function
    at all `Input`).

  - A list of functions. Each function is associated with one cluster.
    These functions are all evaluated at all `Input` values, to provide
    specific hyper-prior mean vectors for each cluster.

- grid_inputs:

  A vector or a data frame, indicating the grid of additional reference
  inputs on which the mean process' hyper-posterior should be evaluated.

- pen_diag:

  A number. A jitter term, added on the diagonal to prevent numerical
  issues when inverting nearly singular matrices.

## Value

A list containing the parameters of the mean processes' hyper-posterior
distribution, namely:

- mean: A list of tibbles containing, for each cluster, the
  hyper-posterior mean parameters evaluated at each `Input`.

- cov: A list of matrices containing, for each cluster, the
  hyper-posterior covariance parameter of the mean process.

- mixture: A tibble, indicating the mixture probabilities in each
  cluster for each individual.

## Examples

``` r
TRUE
#> [1] TRUE
```
