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
#>  1     0   0.97895   0.08123   0.17404  -0.06901
#>  2     1   0.97932   0.07622   0.17471  -0.06786
#>  3     2   0.97988   0.07342   0.17175  -0.07036
#>  4     3   0.98156   0.05804   0.16650  -0.07385
#>  5     4   0.98150   0.06548   0.16639  -0.06841
#>  6     5   0.98178   0.06012   0.16649  -0.06905
#>  7     6   0.98323   0.05813   0.15478  -0.07698
#>  8     7   0.98279   0.06071   0.15564  -0.07880
#>  9     8   0.98356   0.05603   0.15091  -0.08185
#> 10     9   0.98320   0.05944   0.14837  -0.08817
#> # ℹ 91 more rows
#> 
#> attr(,"class")
#> [1] "qts_sample" "list"      
```
