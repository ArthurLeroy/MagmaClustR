# Learning hyper-parameters of a Gaussian Process

Learning hyper-parameters of any new individual/task in `Magma` is
required in the prediction procedure. This function can also be used to
learn hyper-parameters of a simple GP (just let the `hyperpost` argument
set to NULL, and use `prior_mean` instead). When using within `Magma`,
by providing data for the new individual/task, the hyper-posterior mean
and covariance parameters, and initialisation values for the
hyper-parameters, the function computes maximum likelihood estimates of
the hyper-parameters.

## Usage

``` r
train_gp(
  data,
  prior_mean = NULL,
  ini_hp = NULL,
  kern = "SE",
  hyperpost = NULL,
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
  `Input`.

- prior_mean:

  Mean parameter of the GP. This argument can be specified under various
  formats, such as:

  - NULL (default). The hyper-posterior mean would be set to 0
    everywhere.

  - A number. The hyper-posterior mean would be a constant function.

  - A vector of the same length as all the distinct Input values in the
    `data` argument. This vector would be considered as the evaluation
    of the hyper-posterior mean function at the training Inputs.

  - A function. This function is defined as the hyper-posterior mean.

  - A tibble or data frame. Required columns: Input, Output. The Input
    values should include at least the same values as in the `data`
    argument.

- ini_hp:

  A named vector, tibble or data frame of hyper-parameters associated
  with the `kern` of the new individual/task. The columns should be
  named according to the hyper-parameters that are used in `kern`. In
  cases where the model includes a noise term, `ini_hp` should contain
  an additional 'noise' column. If NULL (default), random values are
  used as initialisation. The
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
    separated by whitespaces (e.g. "SE + PERIO"). As the² elements are
    treated sequentially from the left to the right, the product
    operator '\*' shall always be used before the '+' operators (e.g.
    'SE \* LIN + RQ' is valid whereas 'RQ + SE \* LIN' is not).

- hyperpost:

  A list, containing the elements 'mean' and 'cov', the parameters of
  the hyper-posterior distribution of the mean process. Typically, this
  argument should come from a previous learning using
  [`train_magma`](https://arthurleroy.github.io/MagmaClustR/reference/train_magma.md),
  or from the
  [`hyperposterior`](https://arthurleroy.github.io/MagmaClustR/reference/hyperposterior.md)
  function. If `hyperpost` is provided, the likelihood that is maximised
  is the one involved during Magma's prediction step, and the
  `prior_mean` argument is ignored. For classic GP training, leave
  `hyperpost` to NULL.

- pen_diag:

  A number. A jitter term, added on the diagonal to prevent numerical
  issues when inverting nearly singular matrices.

## Value

A tibble, containing the trained hyper-parameters for the kernel of the
new individual/task.

## Examples

``` r
TRUE
#> [1] TRUE
```
