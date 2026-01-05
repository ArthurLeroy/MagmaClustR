# Create inverse of a covariance matrix from a kernel

`kern_to_inv()` creates the inverse of a covariance matrix between input
values (that could be either scalars or vectors) evaluated within a
kernel function, which is characterised by specified hyper-parameters.
This matrix is a finite-dimensional evaluation of the
infinite-dimensional covariance structure of a GP, defined thanks to
this kernel.

## Usage

``` r
kern_to_inv(input, kern, hp, pen_diag = 1e-10, deriv = NULL)
```

## Arguments

- input:

  A vector, matrix, data frame or tibble containing all inputs for one
  individual. If a vector, the elements are used as reference, otherwise
  ,one column should be named 'Input' to indicate that it represents the
  reference (e.g. 'Input' would contain the timestamps in time-series
  applications). The other columns are considered as being covariates.
  If no column is named 'Input', the first one is used by default.

- kern:

  A kernel function. Several popular kernels (see [The Kernel
  Cookbook](https://www.cs.toronto.edu/~duvenaud/cookbook/)) are already
  implemented and can be selected within the following list:

  - "SE": (default value) the Squared Exponential Kernel (also called
    Radial Basis Function or Gaussian kernel),

  - "LIN": the Linear kernel,

  - "PERIO": the Periodic kernel,

  - "RQ": the Rational Quadratic kernel. Compound kernels can be created
    as sums or products of the above kernels. For combining kernels,
    simply provide a formula as a character string where elements are
    separated by whitespaces (e.g. "SE + PERIO"). As the elements are
    treated sequentially from the left to the right, the product
    operator '\*' shall always be used before the '+' operators (e.g.
    'SE \* LIN + RQ' is valid whereas 'RQ + SE \* LIN' is not).

- hp:

  A list, data frame or tibble containing the hyper-parameters used in
  the kernel. The name of the elements (or columns) should correspond
  exactly to those used in the kernel definition.

- pen_diag:

  A jitter term that is added to the covariance matrix to avoid
  numerical issues when inverting, in cases of nearly singular matrices.

- deriv:

  A character, indicating according to which hyper-parameter the
  derivative should be computed. If NULL (default), the function simply
  returns the inverse covariance matrix.

## Value

The inverse of a covariance matrix, which elements are evaluations of
the associated kernel for each pair of reference inputs.

## Examples

``` r
TRUE
#> [1] TRUE
```
