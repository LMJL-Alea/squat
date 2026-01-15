# QTS Transformation to Angular Velocity Vector Time Series

This function projects a quaternion time series into the space of
angular velocity vectors.

## Usage

``` r
qts2avvts(x, spar = 0)
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
with columns `time`, `x`, `y` and `z` containing the angular velocity
vector at each time point.

## Examples

``` r
qts2avvts(vespa64$igp[[1]])
#> # A tibble: 101 × 4
#>     time        x        y         z
#>    <int>    <dbl>    <dbl>     <dbl>
#>  1     0 -0.0109  -0.00444 -0.000474
#>  2     1 -0.0103  -0.00519 -0.000562
#>  3     2 -0.0109  -0.00700 -0.00129 
#>  4     3 -0.00994 -0.00728 -0.00201 
#>  5     4 -0.00832 -0.00690 -0.00255 
#>  6     5 -0.00674 -0.00664 -0.00298 
#>  7     6 -0.00514 -0.00640 -0.00349 
#>  8     7 -0.00383 -0.00642 -0.00405 
#>  9     8 -0.00313 -0.00692 -0.00463 
#> 10     9 -0.00324 -0.00747 -0.00496 
#> # ℹ 91 more rows
```
