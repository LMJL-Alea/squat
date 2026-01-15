# QTS Centering and Standardization

This function operates a centering of the QTS around the geometric mean
of its quaternions. This is effectively achieved by left-multiplying
each quaternion by the inverse of their geometric mean.

## Usage

``` r
centring(x, standardize = FALSE, keep_summary_stats = FALSE)
```

## Arguments

- x:

  An object of class
  [`qts`](https://lmjl-alea.github.io/squat/reference/qts.md).

- standardize:

  A boolean specifying whether to standardize the QTS in addition to
  centering it. Defaults to `FALSE`.

- keep_summary_stats:

  A boolean specifying whether the mean and standard deviation used for
  standardizing the data should be stored in the output object. Defaults
  to `FALSE` in which case only the centered
  [`qts`](https://lmjl-alea.github.io/squat/reference/qts.md) is
  returned.

## Value

If `keep_summary_stats = FALSE`, an object of class
[`qts`](https://lmjl-alea.github.io/squat/reference/qts.md) in which
quaternions have been centered (and possibly standardized) around their
geometric mean. If `keep_summary_stats = TRUE`, a list with three
components:

- `qts`: an object of class
  [`qts`](https://lmjl-alea.github.io/squat/reference/qts.md) in which
  quaternions have been centered (and possibly standardized) around
  their geometric mean;

- `mean`: a numeric vector with the quaternion Fréchet mean;

- `sd`: a numeric value with the quaternion Fréchet standard deviation.

## Examples

``` r
centring(vespa64$igp[[1]])
#> # A tibble: 101 × 5
#>     time         w         x         y         z
#>    <int> <dec:.5!> <dec:.5!> <dec:.5!> <dec:.5!>
#>  1     0   0.99364   0.09383   0.06049   0.01484
#>  2     1   0.99425   0.08868   0.05824   0.01455
#>  3     2   0.99489   0.08343   0.05518   0.01401
#>  4     3   0.99553   0.07811   0.05153   0.01307
#>  5     4   0.99606   0.07365   0.04805   0.01184
#>  6     5   0.99650   0.06989   0.04468   0.01037
#>  7     6   0.99686   0.06694   0.04146   0.00869
#>  8     7   0.99714   0.06478   0.03830   0.00674
#>  9     8   0.99738   0.06307   0.03503   0.00451
#> 10     9   0.99761   0.06156   0.03147   0.00204
#> # ℹ 91 more rows
```
