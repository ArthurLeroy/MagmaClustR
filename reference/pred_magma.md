# Magma prediction

Compute the posterior predictive distribution in Magma. Providing data
of any new individual/task, its trained hyper-parameters and a
previously trained Magma model, the predictive distribution is evaluated
on any arbitrary inputs that are specified through the 'grid_inputs'
argument.

## Usage

``` r
pred_magma(
  data = NULL,
  trained_model = NULL,
  grid_inputs = NULL,
  hp = NULL,
  kern = "SE",
  hyperpost = NULL,
  get_hyperpost = FALSE,
  get_full_cov = FALSE,
  plot = TRUE,
  pen_diag = 1e-10
)
```

## Arguments

- data:

  A tibble or data frame. Required columns: 'Input', 'Output'.
  Additional columns for covariates can be specified. The 'Input' column
  should define the variable that is used as reference for the
  observations (e.g. time for longitudinal data). The 'Output' column
  specifies the observed values (the response variable). The data frame
  can also provide as many covariates as desired, with no constraints on
  the column names. These covariates are additional inputs (explanatory
  variables) of the models that are also observed at each reference
  'Input'. If NULL, the mean process from `trained_model` is returned as
  a generic prediction.

- trained_model:

  A list, containing the information coming from a Magma model,
  previously trained using the
  [`train_magma`](https://arthurleroy.github.io/MagmaClustR/reference/train_magma.md)
  function.

- grid_inputs:

  The grid of inputs (reference Input and covariates) values on which
  the GP should be evaluated. Ideally, this argument should be a tibble
  or a data frame, providing the same columns as `data`, except
  'Output'. Nonetheless, in cases where `data` provides only one 'Input'
  column, the `grid_inputs` argument can be NULL (default) or a vector.
  This vector would be used as reference input for prediction and if
  NULL, a vector of length 500 is defined, ranging between the min and
  max Input values of `data`.

- hp:

  A named vector, tibble or data frame of hyper-parameters associated
  with `kern`. The columns/elements should be named according to the
  hyper-parameters that are used in `kern`. The function
  [`train_gp`](https://arthurleroy.github.io/MagmaClustR/reference/train_gp.md)
  can be used to learn maximum-likelihood estimators of the
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

  A list, containing the elements 'mean' and 'cov', the parameters of
  the hyper-posterior distribution of the mean process. Typically, this
  argument should come from a previous learning using
  [`train_magma`](https://arthurleroy.github.io/MagmaClustR/reference/train_magma.md),
  or a previous prediction with `pred_magma`, with the argument
  `get_hyperpost` set to TRUE. The 'mean' element should be a data frame
  with two columns 'Input' and 'Output'. The 'cov' element should be a
  covariance matrix with colnames and rownames corresponding to the
  'Input' in 'mean'. In all cases, the column 'Input' should contain all
  the values appearing both in the 'Input' column of `data` and in
  `grid_inputs`.

- get_hyperpost:

  A logical value, indicating whether the hyper-posterior distribution
  of the mean process should be returned. This can be useful when
  planning to perform several predictions on the same grid of inputs,
  since recomputation of the hyper-posterior can be prohibitive for high
  dimensional grids.

- get_full_cov:

  A logical value, indicating whether the full posterior covariance
  matrix should be returned.

- plot:

  A logical value, indicating whether a plot of the results is
  automatically displayed.

- pen_diag:

  A number. A jitter term, added on the diagonal to prevent numerical
  issues when inverting nearly singular matrices.

## Value

A tibble, representing Magma predictions as two column 'Mean' and 'Var',
evaluated on the `grid_inputs`. The column 'Input' and additional
covariates columns are associated to each predicted values. If the
`get_full_cov` or `get_hyperpost` arguments are TRUE, the function
returns a list, in which the tibble described above is defined as
'pred_gp' and the full posterior covariance matrix is defined as 'cov',
and the hyper-posterior distribution of the mean process is defined as
'hyperpost'.

## Examples

``` r
TRUE
#> [1] TRUE
```
