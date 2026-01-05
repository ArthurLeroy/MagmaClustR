# Mixture initialisation with kmeans

Provide an initial kmeans allocation of the individuals/tasks in a
dataset into a definite number of clusters, and return the associated
mixture probabilities.

## Usage

``` r
ini_mixture(data, k, name_clust = NULL, nstart = 50)
```

## Arguments

- data:

  A tibble or data frame. Required columns: `ID`, `Input` , `Output`.

- k:

  A number, indicating the number of clusters.

- name_clust:

  A vector of characters. Each element should correspond to the name of
  one cluster.

- nstart:

  A number of restart used in the underlying kmeans algorithm

## Value

A tibble indicating for each `ID` in which cluster it belongs after a
kmeans initialisation.

## Examples

``` r
TRUE
#> [1] TRUE
```
