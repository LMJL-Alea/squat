# QTS Geometric Median

This function computes the pointwise geometric median of a QTS sample.

## Usage

``` r
# S3 method for class 'qts_sample'
median(x, na.rm = FALSE, ...)
```

## Arguments

- x:

  An object of class
  [`qts_sample`](https://lmjl-alea.github.io/squat/reference/qts_sample.md).

- na.rm:

  A logical value indicating whether NA values should be stripped before
  the computation proceeds.

- ...:

  Further arguments passed to or from other methods.

## Value

An object of class
[`qts`](https://lmjl-alea.github.io/squat/reference/qts.md) in which
quaternions are the pointwise geometric median of the input QTS sample.

## Examples

``` r
median(vespa64$igp)
#> # A tibble: 101 × 5
#>     time         w         x         y         z
#>    <int> <dec:.5!> <dec:.5!> <dec:.5!> <dec:.5!>
#>  1     0   0.99394   0.09157   0.05888   0.01513
#>  2     1   0.99439   0.08642   0.05913   0.01490
#>  3     2   0.99482   0.08181   0.05855   0.01465
#>  4     3   0.99525   0.07752   0.05723   0.01414
#>  5     4   0.99565   0.07364   0.05544   0.01338
#>  6     5   0.99603   0.07017   0.05331   0.01236
#>  7     6   0.99638   0.06712   0.05089   0.01105
#>  8     7   0.99671   0.06452   0.04823   0.00938
#>  9     8   0.99701   0.06221   0.04533   0.00742
#> 10     9   0.99729   0.05996   0.04228   0.00530
#> # ℹ 91 more rows
```
