# QTS Transformation to Angular Acceleration Vector Time Series

This function projects a quaternion time series into the space of
angular acceleration vectors.

## Usage

``` r
qts2aavts(x, spar = 0)
```

## Arguments

- x:

  An object of class
  [qts](https://lmjl-alea.github.io/squat/reference/qts.md).

- spar:

  smoothing parameter, typically (but not necessarily) in \\(0,1\]\\.
  When `spar` is specified, the coefficient \\\lambda\\ of the integral
  of the squared second derivative in the fit (penalized log likelihood)
  criterion is a monotone function of `spar`, see the details below.
  Alternatively `lambda` may be specified instead of the *scale free*
  `spar`=\\s\\.

## Value

A time series stored as a
[tibble::tibble](https://tibble.tidyverse.org/reference/tibble.html)
with columns `time`, `x`, `y` and `z` containing the angular
acceleration vector at each time point.

## Examples

``` r
qts2aavts(vespa64$igp[[1]])
#> # A tibble: 101 × 4
#>     time         x         y          z
#>    <int>     <dbl>     <dbl>      <dbl>
#>  1     0  0.00234   0.00106   0.000561 
#>  2     1 -0.00128  -0.00256  -0.000736 
#>  3     2  0.000199 -0.00104  -0.000725 
#>  4     3  0.00168   0.000476 -0.000713 
#>  5     4  0.00156   0.000278 -0.000373 
#>  6     5  0.00159   0.000252 -0.000467 
#>  7     6  0.00162   0.000225 -0.000560 
#>  8     7  0.00100  -0.000258 -0.000572 
#>  9     8  0.000387 -0.000741 -0.000583 
#> 10     9 -0.000597 -0.000377 -0.0000683
#> # ℹ 91 more rows
```
