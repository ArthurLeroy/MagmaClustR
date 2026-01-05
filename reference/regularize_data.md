# Regularise a grid of inputs in a dataset

Modify the original grid of inputs to make it more 'regular' (in the
sense that the interval between each observation is constant, or
corresponds to a specific pattern defined by the user). In particular,
this function can also be used to summarise several data points into
one, at a specific location. In this case, the output values are
averaged according to the 'summarise_fct' argument.

## Usage

``` r
regularize_data(
  data,
  size_grid = 30,
  grid_inputs = NULL,
  summarise_fct = base::mean
)

regularise_data(
  data,
  size_grid = 30,
  grid_inputs = NULL,
  summarise_fct = base::mean
)
```

## Arguments

- data:

  A tibble or data frame. Required columns: `ID`, `Input` `Output`. The
  `ID` column contains the unique names/codes used to identify each
  individual/task (or batch of data). The `Input` column corresponds to
  observed locations (an explanatory variable). The `Output` column
  specifies the associated observed values (the response variable). The
  data frame can also provide as many additional inputs as desired, with
  no constraints on the column names.

- size_grid:

  An integer, which indicates the number of equispaced points each
  column must contain. Each original input value will be collapsed to
  the closest point of the new regular grid, and the associated outputs
  are averaged using the 'summarise_fct' function. This argument is used
  when 'grid_inputs' is left to 'NULL'. Default value is 30.

- grid_inputs:

  A data frame, corresponding to a pre-defined grid of inputs according
  to which we want to regularise a dataset. Column names must be similar
  to those appearing in `data`. If NULL (default), a default grid of
  inputs is defined: for each input column in `data`, a regular sequence
  is created from the min to the max values, with a number of equispaced
  points being equal to the 'size_grid' argument.

- summarise_fct:

  A character string or a function. If several similar inputs are
  associated with different outputs, the user can choose the summarising
  function for the output among the following: min, max, mean, median. A
  custom function can be defined if necessary. Default is "mean".

## Value

A data frame, where input columns have been regularised as desired.

## Examples

``` r
data = tibble::tibble(ID = 1, Input = 0:100, Output = -50:50)

## Define a 1D input grid of 10 points
regularize_data(data, size_grid = 10)
#> # A tibble: 10 × 3
#>       ID Input Output
#>    <dbl> <dbl>  <dbl>
#>  1     1   0    -47.5
#>  2     1  11.1  -39  
#>  3     1  22.2  -28  
#>  4     1  33.3  -17  
#>  5     1  44.4   -5.5
#>  6     1  55.6    6  
#>  7     1  66.7   17  
#>  8     1  77.8   28  
#>  9     1  88.9   39  
#> 10     1 100     47.5

## Define a 1D custom grid
my_grid = tibble::tibble(Input = c(5, 10, 25, 50, 100))
regularize_data(data, grid_inputs = my_grid)
#> # A tibble: 5 × 3
#>      ID Input Output
#>   <dbl> <dbl>  <dbl>
#> 1     1     5  -46.5
#> 2     1    10  -37.5
#> 3     1    25  -22.5
#> 4     1    50    6  
#> 5     1   100   37.5

## Define a 2D input grid of 5x5 points
data_2D = cbind(ID = 1, expand.grid(Input=1:10, Input2=1:10), Output = 1:100)
regularize_data(data_2D, size_grid = 5)
#> # A tibble: 25 × 4
#>       ID Input Input2 Output
#>    <dbl> <dbl>  <dbl>  <dbl>
#>  1     1  0      0       1  
#>  2     1  0      2.25   16  
#>  3     1  0      4.5    36  
#>  4     1  0      6.75   56  
#>  5     1  0      9      81  
#>  6     1  2.25   0       2.5
#>  7     1  2.25   2.25   17.5
#>  8     1  2.25   4.5    37.5
#>  9     1  2.25   6.75   57.5
#> 10     1  2.25   9      82.5
#> # ℹ 15 more rows

## Define a 2D custom input grid
my_grid_2D = MagmaClustR::expand_grid_inputs(c(2, 4, 8), 'Input2' = c(3, 5))
regularize_data(data_2D, grid_inputs = my_grid_2D)
#> # A tibble: 6 × 4
#>      ID Input Input2 Output
#>   <dbl> <dbl>  <dbl>  <dbl>
#> 1     1     2      3   11.5
#> 2     1     2      5   61.5
#> 3     1     4      3   14  
#> 4     1     4      5   64  
#> 5     1     8      3   18  
#> 6     1     8      5   68  
```
