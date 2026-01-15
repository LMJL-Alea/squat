# QTS Transformation to Roll-Pitch-Yaw Time Series

This function converts a quaternion time series into its roll-pitch-yaw
angles representation.

## Usage

``` r
qts2rpyts(x)
```

## Arguments

- x:

  An object of class
  [qts](https://lmjl-alea.github.io/squat/reference/qts.md).

## Value

A time series stored as a
[tibble::tibble](https://tibble.tidyverse.org/reference/tibble.html)
with columns `time`, `roll`, `pitch` and `yaw` containing the roll,
pitch and yaw angles representation of the input quaternions.

## Examples

``` r
qts2rpyts(vespa64$igp[[1]])
#> # A tibble: 101 × 4
#>     time   roll  pitch     yaw
#>    <int>  <dbl>  <dbl>   <dbl>
#>  1     0 0.163  0.137  0.0380 
#>  2     1 0.152  0.133  0.0365 
#>  3     2 0.141  0.127  0.0345 
#>  4     3 0.130  0.120  0.0316 
#>  5     4 0.121  0.113  0.0284 
#>  6     5 0.113  0.107  0.0248 
#>  7     6 0.107  0.101  0.0209 
#>  8     7 0.102  0.0945 0.0165 
#>  9     8 0.0983 0.0882 0.0117 
#> 10     9 0.0950 0.0812 0.00638
#> # ℹ 91 more rows
```
