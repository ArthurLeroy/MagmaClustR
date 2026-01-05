# Squared Exponential Kernel

Squared Exponential Kernel

## Usage

``` r
se_kernel(x, y, hp, deriv = NULL, vectorized = FALSE)
```

## Arguments

- x:

  A vector (or matrix if vectorized = T) of inputs.

- y:

  A vector (or matrix if vectorized = T) of inputs.

- hp:

  A tibble, data frame or named vector, containing the kernel's
  hyperparameters. Required columns: 'se_variance', 'se_lengthscale'.

- deriv:

  A character, indicating according to which hyper-parameter the
  derivative should be computed. If NULL (default), the function simply
  returns the evaluation of the kernel.

- vectorized:

  A logical value, indicating whether the function provides a vectorized
  version for speeded-up calculations. If TRUE, the `x` and `y`
  arguments should be the vector or matrix containing all inputs for
  which the kernel is evaluated on all pairs of elements. If FALSE, the
  `x` and `y` arguments are simply two inputs.

## Value

A scalar, corresponding to the evaluation of the kernel.

## Examples

``` r
TRUE
#> [1] TRUE
```
