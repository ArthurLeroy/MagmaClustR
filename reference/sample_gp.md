# Draw samples from a posterior GP/Magma distribution

Draw samples from a posterior GP/Magma distribution

## Usage

``` r
sample_gp(pred_gp, nb_samples = 50)

sample_magma(pred_gp, nb_samples = 50)
```

## Arguments

- pred_gp:

  A list, typically coming from
  [`pred_magma`](https://arthurleroy.github.io/MagmaClustR/reference/pred_magma.md)
  or
  [`pred_gp`](https://arthurleroy.github.io/MagmaClustR/reference/pred_gp.md)
  functions, with argument 'get_full_cov = TRUE'. Required elements:
  `pred`, `cov`.

- nb_samples:

  A number, indicating the number of samples to be drawn from the
  predictive posterior distribution. For two-dimensional graphs, only
  one sample can be displayed.

## Value

A tibble or data frame, containing the samples generated from a GP
prediction. Format: `Input`, `Sample`, `Output`.

## Examples

``` r
TRUE
#> [1] TRUE
```
