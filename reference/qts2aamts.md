# QTS Transformation to Angular Acceleration Magnitude Time Series

This function projects a quaternion time series into the space of
angular acceleration magnitudes.

## Usage

``` r
qts2aamts(x, spar = 0)
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
with columns `time` and `mag` containing the angular acceleration
magnitude at each time point.

## Examples

``` r
qts2aamts(vespa64$igp[[1]])
#> # A tibble: 101 × 2
#>     time magnitude
#>    <int>     <dbl>
#>  1     0  0.00263 
#>  2     1  0.00296 
#>  3     2  0.00129 
#>  4     3  0.00189 
#>  5     4  0.00163 
#>  6     5  0.00168 
#>  7     6  0.00173 
#>  8     7  0.00118 
#>  9     8  0.00102 
#> 10     9  0.000709
#> # ℹ 91 more rows
```
