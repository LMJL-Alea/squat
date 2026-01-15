# QTS Transformation To Norm Time Series

This function computes a univariate time series representing the norm of
the quaternions.

## Usage

``` r
qts2nts(x, disable_normalization = FALSE)
```

## Arguments

- x:

  An object of class
  [qts](https://lmjl-alea.github.io/squat/reference/qts.md).

- disable_normalization:

  A boolean specifying whether quaternion normalization should be
  disabled. Defaults to `FALSE`.

## Value

A time series stored as a
[tibble::tibble](https://tibble.tidyverse.org/reference/tibble.html)
with columns `time` and `norm` in which `norm` measures the angular
distance between the current quaternion and the identity.

## Examples

``` r
qts2nts(vespa64$igp[[1]])
#> # A tibble: 101 × 2
#>     time  norm
#>    <int> <dbl>
#>  1     0 0.214
#>  2     1 0.203
#>  3     2 0.191
#>  4     3 0.178
#>  5     4 0.167
#>  6     5 0.157
#>  7     6 0.147
#>  8     7 0.140
#>  9     8 0.132
#> 10     9 0.125
#> # ℹ 91 more rows
```
