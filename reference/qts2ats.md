# QTS Transformation To Angle Time Series

This function computes a univariate time series representing the angle
between the first and other attitudes.

## Usage

``` r
qts2ats(x, disable_normalization = FALSE)
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
with columns `time` and `angle` in which `angle` measures the angle
between the current rotation and the first one.

## Examples

``` r
qts2ats(vespa64$igp[[1]])
#> # A tibble: 101 × 2
#>     time  angle
#>    <int>  <dbl>
#>  1     0 0     
#>  2     1 0.0113
#>  3     2 0.0235
#>  4     3 0.0366
#>  5     4 0.0480
#>  6     5 0.0583
#>  7     6 0.0673
#>  8     7 0.0752
#>  9     8 0.0828
#> 10     9 0.0908
#> # ℹ 91 more rows
```
