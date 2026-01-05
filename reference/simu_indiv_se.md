# Simulate a batch of data

Simulate a batch of output data, corresponding to one individual, coming
from a GP with a the Squared Exponential kernel as covariance structure,
and specified hyper-parameters and input.

## Usage

``` r
simu_indiv_se(ID, input, mean, v, l, sigma)
```

## Arguments

- ID:

  An identification code, whether numeric or character.

- input:

  A vector of numbers. The input variable that is used as 'reference'
  for input and outputs.

- mean:

  A vector of numbers. Prior mean values of the GP.

- v:

  A number. The variance hyper-parameter of the SE kernel.

- l:

  A number. The lengthscale hyper-parameter of the SE kernel.

- sigma:

  A number. The noise hyper-parameter.

## Value

A tibble containing a batch of output data along with input and
additional information for a simulated individual.

## Examples

``` r
TRUE
#> [1] TRUE
```
