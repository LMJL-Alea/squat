# Inverse Operator for `qts` Objects

This function implements the pointwise inverse of a quaternion time
series.

## Usage

``` r
inverse_qts(x)
```

## Arguments

- x:

  An object of class
  [`qts`](https://lmjl-alea.github.io/squat/reference/qts.md).

## Value

An object of class
[`qts`](https://lmjl-alea.github.io/squat/reference/qts.md) storing the
inverse of `x`.

## Examples

``` r
inverse_qts(vespa64$igp[[1]])
#> # A tibble: 101 × 5
#>     time         w         x         y         z
#>    <int> <dec:.5!> <dec:.5!> <dec:.5!> <dec:.5!>
#>  1     0   0.99427  -0.07973  -0.06988  -0.01334
#>  2     1   0.99483  -0.07457  -0.06763  -0.01313
#>  3     2   0.99542  -0.06931  -0.06457  -0.01269
#>  4     3   0.99602  -0.06398  -0.06091  -0.01184
#>  5     4   0.99652  -0.05949  -0.05742  -0.01070
#>  6     5   0.99694  -0.05572  -0.05403  -0.00932
#>  7     6   0.99729  -0.05275  -0.05079  -0.00772
#>  8     7   0.99757  -0.05056  -0.04761  -0.00583
#>  9     8   0.99782  -0.04883  -0.04430  -0.00366
#> 10     9   0.99805  -0.04730  -0.04071  -0.00125
#> # ℹ 91 more rows
```
