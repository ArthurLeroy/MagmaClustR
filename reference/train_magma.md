# Training Magma with an EM algorithm

The hyper-parameters and the hyper-posterior distribution involved in
Magma can be learned thanks to an EM algorithm implemented in
`train_magma`. By providing a dataset, the model hypotheses (hyper-prior
mean parameter and covariance kernels) and initialisation values for the
hyper-parameters, the function computes maximum likelihood estimates of
the HPs as well as the mean and covariance parameters of the Gaussian
hyper-posterior distribution of the mean process.

## Usage

``` r
train_magma(
  data,
  prior_mean = NULL,
  ini_hp_0 = NULL,
  ini_hp_i = NULL,
  kern_0 = "SE",
  kern_i = "SE",
  common_hp = TRUE,
  grid_inputs = NULL,
  pen_diag = 1e-10,
  n_iter_max = 25,
  cv_threshold = 0.001,
  fast_approx = FALSE
)
```

## Arguments

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
  observed at each reference `Input`.

- prior_mean:

  Hyper-prior mean parameter (m_0) of the mean GP. This argument can be
  specified under various formats, such as:

  - NULL (default). The hyper-prior mean would be set to 0 everywhere.

  - A number. The hyper-prior mean would be a constant function.

  - A vector of the same length as all the distinct Input values in the
    `data` argument. This vector would be considered as the evaluation
    of the hyper-prior mean function at the training Inputs.

  - A function. This function is defined as the hyper_prior mean.

  - A tibble or data frame. Required columns: Input, Output. The Input
    values should include at least the same values as in the `data`
    argument.

- ini_hp_0:

  A named vector, tibble or data frame of hyper-parameters associated
  with `kern_0`, the mean process' kernel. The columns/elements should
  be named according to the hyper-parameters that are used in `kern_0`.
  If NULL (default), random values are used as initialisation. The
  [`hp`](https://arthurleroy.github.io/MagmaClustR/reference/hp.md)
  function can be used to draw custom hyper-parameters with the correct
  format.

- ini_hp_i:

  A tibble or data frame of hyper-parameters associated with `kern_i`,
  the individual processes' kernel. Required column : `ID`. The `ID`
  column contains the unique names/codes used to identify each
  individual/task. The other columns should be named according to the
  hyper-parameters that are used in `kern_i`. Compared to `ini_hp_0`
  should contain an additional 'noise' column to initialise the noise
  hyper-parameter of the model. If NULL (default), random values are
  used as initialisation. The
  [`hp`](https://arthurleroy.github.io/MagmaClustR/reference/hp.md)
  function can be used to draw custom hyper-parameters with the correct
  format.

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

- kern_i:

  A kernel function, associated with the individual GPs. ("SE", "PERIO"
  and "RQ" are also available here).

- common_hp:

  A logical value, indicating whether the set of hyper-parameters is
  assumed to be common to all individuals.

- grid_inputs:

  A vector, indicating the grid of additional reference inputs on which
  the mean process' hyper-posterior should be evaluated.

- pen_diag:

  A number. A jitter term, added on the diagonal to prevent numerical
  issues when inverting nearly singular matrices.

- n_iter_max:

  A number, indicating the maximum number of iterations of the EM
  algorithm to proceed while not reaching convergence.

- cv_threshold:

  A number, indicating the threshold of the likelihood gain under which
  the EM algorithm will stop. The convergence condition is defined as
  the difference of likelihoods between two consecutive steps, divided
  by the absolute value of the last one ( \\(LL_n - LL_n-1) / \|LL_n\|\\
  ).

- fast_approx:

  A boolean, indicating whether the EM algorithm should stop after only
  one iteration of the E-step. This advanced feature is mainly used to
  provide a faster approximation of the model selection procedure, by
  preventing any optimisation over the hyper-parameters.

## Value

A list, gathering the results of the EM algorithm used for training in
Magma. The elements of the list are:

- hp_0: A tibble of the trained hyper-parameters for the mean process'
  kernel.

- hp_i: A tibble of all the trained hyper-parameters for the individual
  processes' kernels.

- hyperpost: A sub-list gathering the parameters of the mean processes'
  hyper-posterior distributions, namely:

  - mean: A tibble, the hyper-posterior mean parameter (`Output`)
    evaluated at each training reference `Input`.

  - cov: A matrix, the covariance parameter for the hyper-posterior
    distribution of the mean process.

  - pred: A tibble, the predicted mean and variance at `Input` for the
    mean process' hyper-posterior distribution under a format that
    allows the direct visualisation as a GP prediction.

- ini_args: A list containing the initial function arguments and values
  for the hyper-prior mean, the hyper-parameters. In particular, if
  those arguments were set to NULL, `ini_args` allows us to retrieve the
  (randomly chosen) initialisations used during training.

- seq_loglikelihood: A vector, containing the sequence of log-likelihood
  values associated with each iteration.

- converged: A logical value indicated whether the EM algorithm
  converged or not.

- training_time: Total running time of the complete training.

## Details

The user can specify custom kernel functions for the argument `kern_0`
and `kern_i`. The hyper-parameters used in the kernel should have
explicit names, and be contained within the `hp` argument. `hp` should
typically be defined as a named vector or a data frame. Although it is
not mandatory for the `train_magma` function to run, gradients can be
provided within kernel function definition. See for example
[`se_kernel`](https://arthurleroy.github.io/MagmaClustR/reference/se_kernel.md)
to create a custom kernel function displaying an adequate format to be
used in Magma.

## Examples

``` r
TRUE
#> [1] TRUE
```
