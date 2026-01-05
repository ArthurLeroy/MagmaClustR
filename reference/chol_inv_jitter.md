# Inverse a matrix using an adaptive jitter term

Inverse a matrix from its Choleski decomposition. If (nearly-)singular,
increase the order of magnitude of the jitter term added to the diagonal
until the matrix becomes non-singular.

## Usage

``` r
chol_inv_jitter(mat, pen_diag)
```

## Arguments

- mat:

  A matrix, possibly singular.

- pen_diag:

  A number, a jitter term to add on the diagonal.

## Value

A matrix, inverse of `mat` plus an adaptive jitter term added on the
diagonal.

## Examples

``` r
TRUE
#> [1] TRUE
```
