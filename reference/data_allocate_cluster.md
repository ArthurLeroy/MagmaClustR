# Allocate training data into the most probable cluster

Allocate training data into the most probable cluster

## Usage

``` r
data_allocate_cluster(trained_model)
```

## Arguments

- trained_model:

  A list, containing the information coming from a MagmaClust model,
  previously trained using the
  [`train_magmaclust`](https://arthurleroy.github.io/MagmaClustR/reference/train_magmaclust.md)
  function.

## Value

The original dataset used to train the MagmaClust model, with additional
'Cluster' and associated 'Proba' columns, indicating the most probable
cluster for each individual/task at the end of the training procedure.

## Examples

``` r
TRUE
#> [1] TRUE
```
