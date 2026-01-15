# QTS Transformation to Angle-Axis Time Series

This function converts a quaternion time series into its angle-axis
representation.

## Usage

``` r
qts2aats(x)
```

## Arguments

- x:

  An object of class
  [qts](https://lmjl-alea.github.io/squat/reference/qts.md).

## Value

A time series stored as a
[tibble::tibble](https://tibble.tidyverse.org/reference/tibble.html)
with columns `time`, `angle`, `ux`, `uy` and `uz` containing the
angle-axis representation of the input quaternions.

## Examples

``` r
qts2aats(vespa64$igp[[1]])
#> # A tibble: 101 × 5
#>     time angle    ux    uy     uz
#>    <int> <dbl> <dbl> <dbl>  <dbl>
#>  1     0 0.214 0.746 0.654 0.125 
#>  2     1 0.203 0.735 0.666 0.129 
#>  3     2 0.191 0.725 0.676 0.133 
#>  4     3 0.178 0.718 0.683 0.133 
#>  5     4 0.167 0.714 0.689 0.128 
#>  6     5 0.157 0.713 0.691 0.119 
#>  7     6 0.147 0.716 0.690 0.105 
#>  8     7 0.140 0.726 0.683 0.0837
#>  9     8 0.132 0.739 0.671 0.0554
#> 10     9 0.125 0.758 0.652 0.0201
#> # ℹ 91 more rows
```
