# Run a k-means algorithm to initialise clusters' allocation

Run a k-means algorithm to initialise clusters' allocation

## Usage

``` r
ini_kmeans(data, k, nstart = 50, summary = FALSE)
```

## Arguments

- data:

  A tibble containing common Input and associated Output values to
  cluster.

- k:

  A number of clusters assumed for running the kmeans algorithm.

- nstart:

  A number, indicating how many re-starts of kmeans are set.

- summary:

  A boolean, indicating whether we want an outcome summary

## Value

A tibble containing the initial clustering obtained through kmeans.

## Examples

``` r
TRUE
#> [1] TRUE
```
