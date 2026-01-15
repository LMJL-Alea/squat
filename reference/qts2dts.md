# QTS Transformation To Distance Time Series

This function computes a real-valued time series reporting the pointwise
geodesic distance between the two input QTS at each time point.

## Usage

``` r
qts2dts(x, y)
```

## Arguments

- x:

  An object of class
  [qts](https://lmjl-alea.github.io/squat/reference/qts.md).

- y:

  An object of class
  [qts](https://lmjl-alea.github.io/squat/reference/qts.md).

## Value

A time series stored as a
[tibble::tibble](https://tibble.tidyverse.org/reference/tibble.html)
with columns `time` and `distance` in which `distance` measures the
angular distance between the quaternions of both input QTS at a given
time point.

## Details

The function currently expects that the two input QTS are evaluated on
the same time grid.

## Examples

``` r
qts2dts(vespa64$igp[[1]], vespa64$igp[[2]])
#> # A tibble: 101 × 2
#>     time distance
#>    <int>    <dbl>
#>  1     0  0.0120 
#>  2     1  0.0106 
#>  3     2  0.00983
#>  4     3  0.00936
#>  5     4  0.00794
#>  6     5  0.00714
#>  7     6  0.00731
#>  8     7  0.00775
#>  9     8  0.00808
#> 10     9  0.00824
#> # ℹ 91 more rows
```
