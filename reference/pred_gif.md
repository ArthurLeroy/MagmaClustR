# Magma prediction for ploting GIFs

Generate a Magma or classic GP prediction under a format that is
compatible with a further GIF visualisation of the results. For a Magma
prediction, either the `trained_model` or `hyperpost` argument is
required. Otherwise, a classic GP prediction is applied and the prior
mean can be specified through the `mean` argument.

## Usage

``` r
pred_gif(
  data,
  trained_model = NULL,
  grid_inputs = NULL,
  hyperpost = NULL,
  mean = NULL,
  hp = NULL,
  kern = "SE",
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
  'Input'.

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

- hyperpost:

  A list, containing the elements 'mean' and 'cov', the parameters of
  the hyper-posterior distribution of the mean process. Typically, this
  argument should from a previous learning using
  [`train_magma`](https://arthurleroy.github.io/MagmaClustR/reference/train_magma.md),
  or a previous prediction with
  [`pred_magma`](https://arthurleroy.github.io/MagmaClustR/reference/pred_magma.md),
  with the argument `get_hyperpost` set to TRUE. The 'mean' element
  should be a data frame with two columns 'Input' and 'Output'. The
  'cov' element should be a covariance matrix with colnames and rownames
  corresponding to the 'Input' in 'mean'. In all cases, the column
  'Input' should contain all the values appearing both in the 'Input'
  column of `data` and in `grid_inputs`.

- mean:

  Mean parameter of the GP. This argument can be specified under various
  formats, such as:

  - NULL (default). The mean would be set to 0 everywhere.

  - A number. The mean would be a constant function.

  - A function. This function is defined as the mean.

  - A tibble or data frame. Required columns: Input, Output. The Input
    values should include at least the same values as in the `data`
    argument.

- hp:

  A named vector, tibble or data frame of hyper-parameters associated
  with `kern`. The columns/elements should be named according to the
  hyper-parameters that are used in `kern`. The function
  [`train_gp`](https://arthurleroy.github.io/MagmaClustR/reference/train_gp.md)
  can be used to learn maximum-likelihood estimators of the
  hyper-parameters,

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

- pen_diag:

  A number. A jitter term, added on the diagonal to prevent numerical
  issues when inverting nearly singular matrices.

## Value

A tibble, representing Magma or GP predictions as two column 'Mean' and
'Var', evaluated on the `grid_inputs`. The column 'Input' and additional
covariates columns are associated to each predicted values. An
additional 'Index' column is created for the sake of GIF creation using
the function
[`plot_gif`](https://arthurleroy.github.io/MagmaClustR/reference/plot_gif.md)

## Examples

``` r
TRUE
#> [1] TRUE
```
