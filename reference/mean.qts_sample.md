# QTS Geometric Mean

This function computes the pointwise geometric mean of a QTS sample.

## Usage

``` r
# S3 method for class 'qts_sample'
mean(x, ...)
```

## Arguments

- x:

  An object of class
  [`qts_sample`](https://lmjl-alea.github.io/squat/reference/qts_sample.md).

- ...:

  Further arguments passed to or from other methods.

## Value

An object of class
[`qts`](https://lmjl-alea.github.io/squat/reference/qts.md) in which
quaternions are the pointwise geometric mean of the input QTS sample.

## Examples

``` r
mean(vespa64$igp)
#> # A tibble: 101 × 5
#>     time         w         x         y         z
#>    <int> <dec:.5!> <dec:.5!> <dec:.5!> <dec:.5!>
#>  1     0   0.99406   0.09402   0.05288   0.01472
#>  2     1   0.99445   0.08937   0.05352   0.01464
#>  3     2   0.99487   0.08469   0.05348   0.01450
#>  4     3   0.99527   0.08033   0.05280   0.01415
#>  5     4   0.99564   0.07651   0.05164   0.01351
#>  6     5   0.99597   0.07328   0.05011   0.01254
#>  7     6   0.99627   0.07059   0.04825   0.01122
#>  8     7   0.99655   0.06832   0.04607   0.00960
#>  9     8   0.99682   0.06631   0.04356   0.00773
#> 10     9   0.99707   0.06439   0.04079   0.00570
#> # ℹ 91 more rows
```
