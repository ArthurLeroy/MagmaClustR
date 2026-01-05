# Draw samples from a MagmaClust posterior distribution

Draw samples from a MagmaClust posterior distribution

## Usage

``` r
sample_magmaclust(pred_clust, nb_samples = 50)
```

## Arguments

- pred_clust:

  A list, typically coming from
  [`pred_magmaclust`](https://arthurleroy.github.io/MagmaClustR/reference/pred_magmaclust.md),
  with argument get_full_cov = TRUE'. Required elements: `pred`, `cov`,
  `mixture`.

- nb_samples:

  A number, indicating the number of samples to be drawn from the
  predictive posterior distribution. For two-dimensional graphs, only
  one sample can be displayed.

## Value

A tibble or data frame, containing the samples generated from a GP
prediction. Format: `Cluster`, `Proba`, `Input`, `Sample`, `Output`.

## Examples

``` r
TRUE
#> [1] TRUE
```
