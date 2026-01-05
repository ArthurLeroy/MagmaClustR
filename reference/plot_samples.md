# Display realisations from a (mixture of) GP prediction

Display samples drawn from the posterior of a GP, Magma or MagmaClust
prediction. According to the dimension of the inputs, the graph may
represent curves or a heatmap.

## Usage

``` r
plot_samples(
  pred = NULL,
  samples = NULL,
  nb_samples = 50,
  x_input = NULL,
  plot_mean = TRUE,
  alpha_samples = 0.3
)
```

## Arguments

- pred:

  A list, typically coming from
  [`pred_gp`](https://arthurleroy.github.io/MagmaClustR/reference/pred_gp.md),
  [`pred_magma`](https://arthurleroy.github.io/MagmaClustR/reference/pred_magma.md)
  or
  [`pred_magmaclust`](https://arthurleroy.github.io/MagmaClustR/reference/pred_magmaclust.md)
  functions, using the argument 'get_full_cov = TRUE'. Required
  elements: `pred`, `cov`. This argument is needed if `samples` is
  missing.

- samples:

  A tibble or data frame, containing the samples generated from a GP,
  Magma, or MagmaClust prediction. Required columns: `Input`, `Sample`,
  `Output`. This argument is needed if `pred` is missing.

- nb_samples:

  A number, indicating the number of samples to be drawn from the
  predictive posterior distribution. For two-dimensional graphs, only
  one sample can be displayed.

- x_input:

  A vector of character strings, indicating which 'column' should be
  displayed in the case of multidimensional inputs. If NULL(default) the
  Input' column is used for the x-axis. If providing a 2-dimensional
  vector, the corresponding columns are used for the x-axis and the
  y-axis.

- plot_mean:

  A logical value, indicating whether the mean prediction should be
  displayed on the graph.

- alpha_samples:

  A number, controlling transparency of the sample curves.

## Value

Graph of samples drawn from a posterior distribution of a GP, Magma, or
MagmaClust prediction.

## Examples

``` r
TRUE
#> [1] TRUE
```
