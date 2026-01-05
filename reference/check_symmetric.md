# Round a matrix to make if symmetric

If a matrix is non-symmetric due to numerical errors, round with a
decreasing number of digits until the matrix becomes symmetric.

## Usage

``` r
check_symmetric(mat, digits = 10)
```

## Arguments

- mat:

  A matrix, possibly non-symmetric.

- digits:

  A number, the starting number of digits to round from if `mat` is not
  symmetric

## Value

A matrix, rounded approximation of `mat` that is symmetric.

## Examples

``` r
TRUE
#> [1] TRUE
```
