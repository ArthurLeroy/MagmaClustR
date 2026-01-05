# Compute the Multivariate Gaussian likelihood

Modification of the function `dmvnorm()` from the package `mvtnorm`,
providing an implementation of the Multivariate Gaussian likelihood.
This version uses inverse of the covariance function as argument instead
of the traditional covariance.

## Usage

``` r
dmnorm(x, mu, inv_Sigma, log = FALSE)
```

## Arguments

- x:

  A vector, containing values the likelihood is evaluated on.

- mu:

  A vector or matrix, specifying the mean parameter.

- inv_Sigma:

  A matrix, specifying the inverse of covariance parameter.

- log:

  A logical value, indicating whether we return the log-likelihood.

## Value

A number, corresponding to the Multivariate Gaussian log-likelihood.

## Examples

``` r
TRUE
#> [1] TRUE
```
