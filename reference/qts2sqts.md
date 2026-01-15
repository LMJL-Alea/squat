# QTS Transformation to Smoothed Quaternion Time Series

This function smooths a given QTS using a b-spline functional
representation.

## Usage

``` r
qts2sqts(x, spar = 0)
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

An object of class `qts` storing the smoothed QTS.

## Examples

``` r
qts2sqts(vespa64$igp[[1]], spar = 0.3)
#> # A tibble: 101 × 5
#>     time         w         x         y         z
#>    <int> <dec:.5!> <dec:.5!> <dec:.5!> <dec:.5!>
#>  1     0   0.99044   0.07964   0.07964   0.07964
#>  2     1   0.99167   0.07437   0.07437   0.07437
#>  3     2   0.99281   0.06910   0.06910   0.06910
#>  4     3   0.99382   0.06406   0.06406   0.06406
#>  5     4   0.99467   0.05955   0.05955   0.05955
#>  6     5   0.99532   0.05579   0.05579   0.05579
#>  7     6   0.99581   0.05283   0.05283   0.05283
#>  8     7   0.99615   0.05059   0.05059   0.05059
#>  9     8   0.99642   0.04881   0.04881   0.04881
#> 10     9   0.99665   0.04720   0.04720   0.04720
#> # ℹ 91 more rows
```
