# Select the optimal number of clusters

In MagmaClust, as for any clustering method, the number K of clusters
has to be provided as an hypothesis of the model. This function
implements a model selection procedure, by maximising a variational BIC
criterion, computed for different values of K. A heuristic for a fast
approximation of the procedure is proposed as well, although the
corresponding models would not be properly trained.

## Usage

``` r
select_nb_cluster(
  data,
  fast_approx = TRUE,
  grid_nb_cluster = 1:10,
  ini_hp_k = NULL,
  ini_hp_i = NULL,
  kern_k = "SE",
  kern_i = "SE",
  plot = TRUE,
  ...
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

- fast_approx:

  A boolean, indicating whether a fast approximation should be used for
  selecting the number of clusters. If TRUE, each Magma or MagmaClust
  model will perform only one E-step of the training, using the same
  fixed values for the hyper-parameters (`ini_hp_k` and `ini_hp_i`, or
  random values if not provided) in all models. The resulting models
  should not be considered as trained, but this approach provides an
  convenient heuristic to avoid a cumbersome model selection procedure.

- grid_nb_cluster:

  A vector of integer, corresponding to grid of values that will be
  tested for the number of clusters.

- ini_hp_k:

  A tibble or data frame of hyper-parameters associated with `kern_k`.
  The [`hp`](https://arthurleroy.github.io/MagmaClustR/reference/hp.md)
  function can be used to draw custom hyper-parameters with the correct
  format.

- ini_hp_i:

  A tibble or data frame of hyper-parameters associated with `kern_i`.
  The [`hp`](https://arthurleroy.github.io/MagmaClustR/reference/hp.md)
  function can be used to draw custom hyper-parameters with the correct
  format.db

- kern_k:

  A kernel function associated to the mean processes.

- kern_i:

  A kernel function associated to the individuals/tasks.

- plot:

  A boolean indicating whether the plot of V-BIC values for all numbers
  of clusters should displayed.

- ...:

  Any additional argument that could be passed to
  [`train_magmaclust`](https://arthurleroy.github.io/MagmaClustR/reference/train_magmaclust.md).

## Value

A list, containing the results of model selection procedure for
selecting the optimal number of clusters thanks to a V-BIC criterion
maximisation. The elements of the list are:

- best_k: An integer, indicating the resulting optimal number of
  clusters

- seq_vbic: A vector, corresponding to the sequence of the V-BIC values
  associated with the models trained for each provided cluster's number
  in `grid_nb_cluster`.

- trained_models: A list, named by associated number of clusters, of
  Magma or MagmaClust models that have been trained (or approximated if
  `fast_approx` = T) during the model selection procedure.

## Examples

``` r
TRUE
#> [1] TRUE
```
