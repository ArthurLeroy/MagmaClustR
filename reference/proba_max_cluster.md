# Indicates the most probable cluster

Indicates the most probable cluster

## Usage

``` r
proba_max_cluster(mixture)
```

## Arguments

- mixture:

  A tibble or data frame containing mixture probabilities.

## Value

A tibble, retaining only the most probable cluster. The column `Cluster`
indicates the the cluster's name whereas `Proba` refers to its
associated probability. If `ID` is initially a column of `mixture`
(optional), the function returns the most probable cluster for all the
different `ID` values.

## Examples

``` r
TRUE
#> [1] TRUE
```
