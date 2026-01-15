# Operator `*` for `qts` Objects

This function implements the pointwise quaternion Hamilton
multiplication between two quaternion time series.

## Usage

``` r
# S3 method for class 'qts'
x * rhs
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
multiplication of the two inputs.

## Examples

``` r
vespa64$igp[[1]] * vespa64$igp[[2]]
#> # A tibble: 101 × 5
#>     time         w         x         y         z
#>    <int> <dec:.5!> <dec:.5!> <dec:.5!> <dec:.5!>
#>  1     0   0.97621   0.16129   0.14303   0.02323
#>  2     1   0.97856   0.15079   0.13832   0.02344
#>  3     2   0.98096   0.14017   0.13241   0.02327
#>  4     3   0.98339   0.12929   0.12539   0.02235
#>  5     4   0.98557   0.11930   0.11828   0.02086
#>  6     5   0.98744   0.11058   0.11124   0.01875
#>  7     6   0.98894   0.10387   0.10465   0.01602
#>  8     7   0.99017   0.09883   0.09812   0.01261
#>  9     8   0.99126   0.09487   0.09128   0.00851
#> 10     9   0.99225   0.09158   0.08394   0.00386
#> # ℹ 91 more rows
```
