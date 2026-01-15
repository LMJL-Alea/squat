# QTS Transformation to Angular Velocity Magnitude Time Series

This function projects a quaternion time series into the space of
angular velocity magnitudes.

## Usage

``` r
qts2avmts(x, spar = 0)
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
with columns `time` and `magnitude` containing the angular velocity
magnitude at each time point.

## Examples

``` r
qts2avmts(vespa64$igp[[1]])
#> # A tibble: 101 × 2
#>     time magnitude
#>    <int>     <dbl>
#>  1     0   0.0117 
#>  2     1   0.0116 
#>  3     2   0.0130 
#>  4     3   0.0125 
#>  5     4   0.0111 
#>  6     5   0.00992
#>  7     6   0.00892
#>  8     7   0.00850
#>  9     8   0.00889
#> 10     9   0.00953
#> # ℹ 91 more rows
```
