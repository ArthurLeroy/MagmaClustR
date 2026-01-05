# Compute the hyper-posterior distribution in Magma

Compute the parameters of the hyper-posterior Gaussian distribution of
the mean process in Magma (similarly to the expectation step of the EM
algorithm used for learning). This hyper-posterior distribution,
evaluated on a grid of inputs provided through the `grid_inputs`
argument, is a key component for making prediction in Magma, and is
required in the function
[`pred_magma`](https://arthurleroy.github.io/MagmaClustR/reference/pred_magma.md).

## Usage

``` r
hyperposterior(
  trained_model = NULL,
  data = NULL,
  hp_0 = NULL,
  hp_i = NULL,
  kern_0 = NULL,
  kern_i = NULL,
  prior_mean = NULL,
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
  `hp_0`, `hp_i`, `kern_0`, and `kern_i` are all required.

- data:

  A tibble or data frame. Required columns: 'Input', 'Output'.
  Additional columns for covariates can be specified. The 'Input' column
  should define the variable that is used as reference for the
  observations (e.g. time for longitudinal data). The 'Output' column
  specifies the observed values (the response variable). The data frame
  can also provide as many covariates as desired, with no constraints on
  the column names. These covariates are additional inputs (explanatory
  variables) of the models that are also observed at each reference
  'Input'. Recovered from `trained_model` if not provided.

- hp_0:

  A named vector, tibble or data frame of hyper-parameters associated
  with `kern_0`. Recovered from `trained_model` if not provided.

- hp_i:

  A tibble or data frame of hyper-parameters associated with `kern_i`.
  Recovered from `trained_model` if not provided.

- kern_0:

  A kernel function, associated with the mean GP. Several popular
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

  A kernel function, associated with the individual GPs. ("SE", "PERIO"
  and "RQ" are aso available here). Recovered from `trained_model` if
  not provided.

- prior_mean:

  Hyper-prior mean parameter of the mean GP. This argument, can be
  specified under various formats, such as:

  - NULL (default). The hyper-prior mean would be set to 0 everywhere.

  - A number. The hyper-prior mean would be a constant function.

  - A vector of the same length as all the distinct Input values in the
    `data` argument. This vector would be considered as the evaluation
    of the hyper-prior mean function at the training Inputs.

  - A function. This function is defined as the hyper-prior mean.

  - A tibble or data frame. Required columns: Input, Output. The Input
    values should include at least the same values as in the `data`
    argument.

- grid_inputs:

  A vector or a data frame, indicating the grid of additional reference
  inputs on which the mean process' hyper-posterior should be evaluated.

- pen_diag:

  A number. A jitter term, added on the diagonal to prevent numerical
  issues when inverting nearly singular matrices.

## Value

A list gathering the parameters of the mean processes' hyper-posterior
distributions, namely:

- mean: A tibble, the hyper-posterior mean parameter evaluated at each
  training `Input`.

- cov: A matrix, the covariance parameter for the hyper-posterior
  distribution of the mean process.

- pred: A tibble, the predicted mean and variance at `Input` for the
  mean process' hyper-posterior distribution under a format that allows
  the direct visualisation as a GP prediction.

## Examples

``` r
TRUE
#> [1] TRUE
```
