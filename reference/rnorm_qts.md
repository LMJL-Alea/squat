# QTS Random Sampling

This function adds uncorrelated Gaussian noise to the logarithm QTS
using an exponential covariance function.

## Usage

``` r
rnorm_qts(n, mean_qts, alpha = 0.01, beta = 0.001)
```

## Arguments

- n:

  An integer specifying how many QTS should be generated.

- mean_qts:

  An object of class
  [`qts`](https://lmjl-alea.github.io/squat/reference/qts.md) specifying
  the mean QTS.

- alpha:

  A positive scalar specifying the variance of each component of the
  log-QTS. Defaults to `0.01`.

- beta:

  A positive scalar specifying the exponential weight. Defaults to
  `0.001`.

## Value

A list of `n` objects of class
[`qts`](https://lmjl-alea.github.io/squat/reference/qts.md) with added
noise as specified by parameters `alpha` and `beta`.

## Details

See
[`exp_cov_function`](https://astamm.github.io/roahd/reference/exp_cov_function.html)
for details about the roles of `alpha` and `beta` in the definition of
the covariance operator.

## Examples

``` r
rnorm_qts(1, vespa64$igp[[1]])
#> [[1]]
#> # A tibble: 101 × 5
#>     time         w         x         y         z
#>    <int> <dec:.5!> <dec:.5!> <dec:.5!> <dec:.5!>
#>  1     0   0.91688   0.33844   0.14417   0.15493
#>  2     1   0.91902   0.33181   0.14549   0.15535
#>  3     2   0.92323   0.32782   0.13435   0.14876
#>  4     3   0.92359   0.32706   0.13479   0.14777
#>  5     4   0.92708   0.31990   0.13011   0.14577
#>  6     5   0.92788   0.32089   0.12790   0.14037
#>  7     6   0.93216   0.31711   0.11702   0.12969
#>  8     7   0.93450   0.30981   0.11048   0.13611
#>  9     8   0.93823   0.30398   0.10235   0.12978
#> 10     9   0.93792   0.30684   0.10035   0.12684
#> # ℹ 91 more rows
#> 
#> attr(,"class")
#> [1] "qts_sample" "list"      
```
