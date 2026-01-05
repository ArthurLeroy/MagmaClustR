# Compute an inverse covariance matrix for multiple individuals

Compute the inverse covariance matrices associated with all individuals
in the database, taking into account their specific inputs and
hyper-parameters.

## Usage

``` r
list_kern_to_inv(db, kern, hp, pen_diag, deriv = NULL)
```

## Arguments

- db:

  A tibble or data frame of input data. Required column: 'ID'. Suggested
  column: 'Input' (for indicating the reference input).

- kern:

  A kernel function.

- hp:

  A tibble or data frame, containing the hyper-parameters associated
  with each individual.

- pen_diag:

  A number. A jitter term, added on the diagonal to prevent numerical
  issues when inverting nearly singular matrices.

- deriv:

  A character, indicating according to which hyper-parameter the
  derivative should be computed. If NULL (default), the function simply
  returns the list of covariance matrices.

## Value

A named list containing all of the inverse covariance matrices.

## Examples

``` r
TRUE
#> [1] TRUE
```
