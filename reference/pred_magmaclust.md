# MagmaClust prediction

Compute the posterior predictive distribution in MagmaClust. Providing
data from any new individual/task, its trained hyper-parameters and a
previously trained MagmaClust model, the multi-task posterior
distribution is evaluated on any arbitrary inputs that are specified
through the 'grid_inputs' argument. Due to the nature of the model, the
prediction is defined as a mixture of Gaussian distributions. Therefore
the present function computes the parameters of the predictive
distribution associated with each cluster, as well as the posterior
mixture probabilities for this new individual/task.

## Usage

``` r
pred_magmaclust(
  data = NULL,
  trained_model = NULL,
  grid_inputs = NULL,
  mixture = NULL,
  hp = NULL,
  kern = "SE",
  hyperpost = NULL,
  prop_mixture = NULL,
  get_hyperpost = FALSE,
  get_full_cov = TRUE,
  plot = TRUE,
  pen_diag = 1e-10
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
  'Input'. If NULL, the mixture of mean processes from `trained_model`
  is returned as a generic prediction.

- trained_model:

  A list, containing the information coming from a MagmaClust model,
  previously trained using the
  [`train_magmaclust`](https://arthurleroy.github.io/MagmaClustR/reference/train_magmaclust.md)
  function. If `trained_model` is set to NULL, the `hyperpost` and
  `prop_mixture` arguments are mandatory to perform required
  re-computations for the prediction to succeed.

- grid_inputs:

  The grid of inputs (reference Input and covariates) values on which
  the GP should be evaluated. Ideally, this argument should be a tibble
  or a data frame, providing the same columns as `data`, except
  'Output'. Nonetheless, in cases where `data` provides only one 'Input'
  column, the `grid_inputs` argument can be NULL (default) or a vector.
  This vector would be used as reference input for prediction and if
  NULL, a vector of length 500 is defined, ranging between the min and
  max Input values of `data`.

- mixture:

  A tibble or data frame, indicating the mixture probabilities of each
  cluster for the new individual/task. If NULL, the
  [`train_gp_clust`](https://arthurleroy.github.io/MagmaClustR/reference/train_gp_clust.md)
  function is used to compute these posterior probabilities according to
  `data`.

- hp:

  A named vector, tibble or data frame of hyper-parameters associated
  with `kern`. The columns/elements should be named according to the
  hyper-parameters that are used in `kern`. The
  [`train_gp_clust`](https://arthurleroy.github.io/MagmaClustR/reference/train_gp_clust.md)
  function can be used to learn maximum-likelihood estimators of the
  hyper-parameters.

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
  or a previous prediction with `pred_magmaclust`, with the argument
  `get_hyperpost` set to TRUE.

- prop_mixture:

  A tibble or a named vector of the mixture proportions. Each name of
  column or element should refer to a cluster. The value associated with
  each cluster is a number between 0 and 1. If both `mixture` and
  `trained_model` are set to NULL, this argument allows to recompute
  mixture probabilities, thanks to the `hyperpost` argument and the
  [`train_gp_clust`](https://arthurleroy.github.io/MagmaClustR/reference/train_gp_clust.md)
  function.

- get_hyperpost:

  A logical value, indicating whether the hyper-posterior distributions
  of the mean processes should be returned. This can be useful when
  planning to perform several predictions on the same grid of inputs,
  since recomputation of the hyper-posterior can be prohibitive for high
  dimensional grids.

- get_full_cov:

  A logical value, indicating whether the full posterior covariance
  matrices should be returned.

- plot:

  A logical value, indicating whether a plot of the results is
  automatically displayed.

- pen_diag:

  A number. A jitter term, added on the diagonal to prevent numerical
  issues when inverting nearly singular matrices.

## Value

A list of GP prediction results composed of:

- pred: As sub-list containing, for each cluster:

  - pred_gp: A tibble, representing the GP predictions as two column
    `Mean` and `Var`, evaluated on the `grid_inputs`. The column `Input`
    and additional covariates columns are associated with each predicted
    values.

  - proba: A number, the posterior probability associated with this
    cluster.

  - cov (if `get_full_cov` = TRUE): A matrix, the full posterior
    covariance matrix associated with this cluster.

- mixture: A tibble, indicating the mixture probabilities of each
  cluster for the predicted individual/task.

- hyperpost (if `get_hyperpost` = TRUE): A list, containing the
  hyper-posterior distributions information useful for visualisation
  purposes.

## Examples

``` r
TRUE
#> [1] TRUE
```
