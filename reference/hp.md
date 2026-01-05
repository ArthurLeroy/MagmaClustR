# Generate random hyper-parameters

Generate a set of random hyper-parameters, specific to the chosen type
of kernel, under the format that is used in Magma.

## Usage

``` r
hp(
  kern = "SE",
  list_ID = NULL,
  list_hp = NULL,
  noise = FALSE,
  common_hp = FALSE
)
```

## Arguments

- kern:

  A function, or a character string indicating the chosen type of kernel
  among:

  - "SE": the Squared Exponential kernel,

  - "LIN": the Linear kernel,

  - "PERIO": the Periodic kernel,

  - "RQ": the Rational Quadratic kernel. Compound kernels can be created
    as sums or products of the above kernels. For combining kernels,
    simply provide a formula as a character string where elements are
    separated by whitespaces (e.g. "SE + PERIO"). As the elements are
    treated sequentially from the left to the right, the product
    operator '\*' shall always be used before the '+' operators (e.g.
    'SE \* LIN + RQ' is valid whereas 'RQ + SE \* LIN' is not).

  In case of a custom kernel function, the argument `list_hp` has to be
  provided as well, for designing a tibble with the correct names of
  hyper-parameters.

- list_ID:

  A vector, associating an `ID` value with each individual for whom
  hyper-parameters are generated. If NULL (default) only one set of
  hyper-parameters is return without the `ID` column.

- list_hp:

  A vector of characters, providing the name of each hyper-parameter, in
  case where `kern` is a custom kernel function.

- noise:

  A logical value, indicating whether a 'noise' hyper-parameter should
  be included.

- common_hp:

  A logical value, indicating whether the set of hyper-parameters is
  assumed to be common to all individuals.

## Value

A tibble, providing a set of random hyper-parameters associated with the
kernel specified through the argument `kern`.

## Examples

``` r
TRUE
#> [1] TRUE
```
