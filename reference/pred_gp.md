# Gaussian Process prediction

Compute the posterior distribution of a standard GP, using the formalism
of Magma. By providing observed data, the prior mean and covariance
matrix (by defining a kernel and its associated hyper-parameters), the
mean and covariance parameters of the posterior distribution are
computed on the grid of inputs that has been specified. This predictive
distribution can be evaluated on any arbitrary inputs since a GP is an
infinite-dimensional object.

## Usage

``` r
pred_gp(
  data = NULL,
  grid_inputs = NULL,
  mean = NULL,
  hp = NULL,
  kern = "SE",
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
  'Input'. If NULL, the prior GP is returned.

- grid_inputs:

  The grid of inputs (reference Input and covariates) values on which
  the GP should be evaluated. Ideally, this argument should be a tibble
  or a data frame, providing the same columns as `data`, except
  'Output'. Nonetheless, in cases where `data` provides only one 'Input'
  column, the `grid_inputs` argument can be NULL (default) or a vector.
  This vector would be used as reference input for prediction and if
  NULL, a vector of length 500 is defined, ranging between the min and
  max Input values of `data`.

- mean:

  Mean parameter of the GP. This argument can be specified under various
  formats, such as:

  - NULL (default). The mean would be set to 0 everywhere.

  - A number. The mean would be a constant function.

  - A tibble or data frame. Required columns: Input, Output. The Input
    values should include at least the same values as in the `data`
    argument.

- hp:

  A named vector, tibble or data frame of hyper-parameters associated
  with `kern`. The columns/elements should be named according to the
  hyper-parameters that are used in `kern`. If NULL (default), the
  function
  [`train_gp`](https://arthurleroy.github.io/MagmaClustR/reference/train_gp.md)
  is called with random initial values for learning maximum-likelihood
  estimators of the hyper-parameters associated with `kern`.

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

A tibble, representing the GP predictions as two column 'Mean' and
'Var', evaluated on the `grid_inputs`. The column 'Input' and additional
covariates columns are associated to each predicted values. If the
`get_full_cov` argument is TRUE, the function returns a list, in which
the tibble described above is defined as 'pred' and the full posterior
covariance matrix is defined as 'cov'.

## Examples

``` r
TRUE
#> [1] TRUE
```
