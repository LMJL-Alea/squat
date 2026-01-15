# Operator `-` for `qts` Objects

This function implements the pointwise subtraction between two
quaternion time series.

## Usage

``` r
# S3 method for class 'qts'
x - rhs
```

## Arguments

- x:

  An object of class
  [`qts`](https://lmjl-alea.github.io/squat/reference/qts.md).

- rhs:

  Either an object of class
  [`qts`](https://lmjl-alea.github.io/squat/reference/qts.md) or a
  numeric value.

## Value

An object of class
[`qts`](https://lmjl-alea.github.io/squat/reference/qts.md) storing the
subtraction of the two inputs.

## Examples

``` r
vespa64$igp[[1]] - vespa64$igp[[2]]
#> # A tibble: 101 × 5
#>     time         w         x         y         z
#>    <int> <dec:.5!> <dec:.5!> <dec:.5!> <dec:.5!>
#>  1     0   0.00049  -0.00309  -0.00382   0.00340
#>  2     1   0.00042  -0.00270  -0.00356   0.00277
#>  3     2   0.00039  -0.00240  -0.00372   0.00209
#>  4     3   0.00037  -0.00201  -0.00398   0.00137
#>  5     4   0.00027  -0.00083  -0.00382   0.00065
#>  6     5   0.00017   0.00047  -0.00353   0.00006
#>  7     6   0.00011   0.00133  -0.00338  -0.00039
#>  8     7   0.00006   0.00207  -0.00319  -0.00072
#>  9     8   0.00001   0.00262  -0.00293  -0.00095
#> 10     9  -0.00001   0.00288  -0.00274  -0.00111
#> # ℹ 91 more rows
```
