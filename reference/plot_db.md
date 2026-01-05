# Plot smoothed curves of raw data

Display raw data under the Magma format as points and lines.

## Usage

``` r
plot_db(data, legend = FALSE)
```

## Arguments

- data:

  A data frame or tibble with format : ID, Input, Output.

- legend:

  A boolean indicating whether the legend should be displayed.

## Value

Graph of the raw data.

## Examples

``` r
## Generate a synthetic dataset
db = simu_db()
## Plot the raw data
plot_db(db)
```
