# Operator `+` for `qts` Objects

This function implements the pointwise addition between two quaternion
time series.

## Usage

``` r
# S3 method for class 'qts'
x + rhs
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
addition of the two inputs.

## Examples

``` r
vespa64$igp[[1]] + vespa64$igp[[2]]
#> # A tibble: 101 × 5
#>     time         w         x         y         z
#>    <int> <dec:.5!> <dec:.5!> <dec:.5!> <dec:.5!>
#>  1     0   1.98806   0.16255   0.14357   0.02328
#>  2     1   1.98924   0.15184   0.13883   0.02348
#>  3     2   1.99045   0.14102   0.13286   0.02328
#>  4     3   1.99167   0.12996   0.12580   0.02232
#>  5     4   1.99277   0.11982   0.11866   0.02076
#>  6     5   1.99371   0.11097   0.11159   0.01859
#>  7     6   1.99446   0.10417   0.10497   0.01582
#>  8     7   1.99508   0.09906   0.09841   0.01238
#>  9     8   1.99562   0.09505   0.09154   0.00827
#> 10     9   1.99611   0.09172   0.08416   0.00362
#> # ℹ 91 more rows
```
