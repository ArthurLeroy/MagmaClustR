# Training MagmaClust with a Variational EM algorithm

The hyper-parameters and the hyper-posterior distributions involved in
MagmaClust can be learned thanks to a VEM algorithm implemented in
`train_magmaclust`. By providing a dataset, the model hypotheses
(hyper-prior mean parameters, covariance kernels and number of clusters)
and initialisation values for the hyper-parameters, the function
computes maximum likelihood estimates of the HPs as well as the mean and
covariance parameters of the Gaussian hyper-posterior distributions of
the mean processes.

## Usage

``` r
train_magmaclust(
  data,
  nb_cluster = NULL,
  prior_mean_k = NULL,
  ini_hp_k = NULL,
  ini_hp_i = NULL,
  kern_k = "SE",
  kern_i = "SE",
  ini_mixture = NULL,
  common_hp_k = TRUE,
  common_hp_i = TRUE,
  grid_inputs = NULL,
  pen_diag = 1e-10,
  n_iter_max = 25,
  cv_threshold = 0.001,
  fast_approx = FALSE
)
```

## Arguments

- data:

  A tibble or data frame. Columns required: `ID`, `Input` , `Output`.
  Additional columns for covariates can be specified. The `ID` column
  contains the unique names/codes used to identify each individual/task
  (or batch of data). The `Input` column should define the variable that
  is used as reference for the observations (e.g. time for longitudinal
  data). The `Output` column specifies the observed values (the response
  variable). The data frame can also provide as many covariates as
  desired, with no constraints on the column names. These covariates are
  additional inputs (explanatory variables) of the models that are also
  observed at each reference `Input`.

- nb_cluster:

  A number, indicating the number of clusters of individuals/tasks that
  are assumed to exist among the dataset.

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

- ini_hp_k:

  A tibble or data frame of hyper-parameters associated with `kern_k`,
  the mean process' kernel. Required column : `ID`. The `ID` column
  contains the unique names/codes used to identify each cluster. The
  other columns should be named according to the hyper-parameters that
  are used in `kern_k`. The
  [`hp`](https://arthurleroy.github.io/MagmaClustR/reference/hp.md)
  function can be used to draw custom hyper-parameters with the correct
  format.

- ini_hp_i:

  A tibble or data frame of hyper-parameters associated with `kern_i`,
  the individual processes' kernel. Required column : `ID`. The `ID`
  column contains the unique names/codes used to identify each
  individual/task. The other columns should be named according to the
  hyper-parameters that are used in `kern_i`. The
  [`hp`](https://arthurleroy.github.io/MagmaClustR/reference/hp.md)
  function can be used to draw custom hyper-parameters with the correct
  format.

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

- kern_i:

  A kernel function, associated with the individual GPs. (See details
  above in `kern_k`).

- ini_mixture:

  Initial values of the probability to belong to each cluster for each
  individual
  ([`ini_mixture`](https://arthurleroy.github.io/MagmaClustR/reference/ini_mixture.md)
  can be used for a k-means initialisation. Used by default if NULL).

- common_hp_k:

  A boolean indicating whether hyper-parameters are common among the
  mean GPs.

- common_hp_i:

  A boolean indicating whether hyper-parameters are common among the
  individual GPs.

- grid_inputs:

  A vector, indicating the grid of additional reference inputs on which
  the mean processes' hyper-posteriors should be evaluated.

- pen_diag:

  A number. A jitter term, added on the diagonal to prevent numerical
  issues when inverting nearly singular matrices.

- n_iter_max:

  A number, indicating the maximum number of iterations of the VEM
  algorithm to proceed while not reaching convergence.

- cv_threshold:

  A number, indicating the threshold of the likelihood gain under which
  the VEM algorithm will stop. The convergence condition is defined as
  the difference of elbo between two consecutive steps, divided by the
  absolute value of the last one ( \\(ELBO_n - ELBO\_{n-1}) / \|ELBO_n\|
  \\ ).

- fast_approx:

  A boolean, indicating whether the VEM algorithm should stop after only
  one iteration of the VE-step. This advanced feature is mainly used to
  provide a faster approximation of the model selection procedure, by
  preventing any optimisation over the hyper-parameters.

## Value

A list, containing the results of the VEM algorithm used in the training
step of MagmaClust. The elements of the list are:

- hp_k: A tibble containing the trained hyper-parameters for the mean
  process' kernel and the mixture proportions for each cluster.

- hp_i: A tibble containing the trained hyper-parameters for the
  individual processes' kernels.

- hyperpost: A sub-list containing the parameters of the mean processes'
  hyper-posterior distribution, namely:

  - mean: A list of tibbles containing, for each cluster, the
    hyper-posterior mean parameters evaluated at each `Input`.

  - cov: A list of matrices containing, for each cluster, the
    hyper-posterior covariance parameter of the mean process.

  - mixture: A tibble, indicating the mixture probabilities in each
    cluster for each individual.

- ini_args: A list containing the initial function arguments and values
  for the hyper-prior means, the hyper-parameters. In particular, if
  those arguments were set to NULL, `ini_args` allows us to retrieve the
  (randomly chosen) initialisations used during training.

- seq_elbo: A vector, containing the sequence of ELBO values associated
  with each iteration.

- converged: A logical value indicated whether the algorithm converged.

- training_time: Total running time of the complete training.

## Details

The user can specify custom kernel functions for the argument `kern_k`
and `kern_i`. The hyper-parameters used in the kernel should have
explicit names, and be contained within the `hp` argument. `hp` should
typically be defined as a named vector or a data frame. Although it is
not mandatory for the `train_magmaclust` function to run, gradients be
can provided within kernel function definition. See for example
[`se_kernel`](https://arthurleroy.github.io/MagmaClustR/reference/se_kernel.md)
to create a custom kernel function displaying an adequate format to be
used in MagmaClust.

## Examples

``` r
TRUE
#> [1] TRUE
```
