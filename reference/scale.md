# QTS Sample Centering and Standardization

QTS Sample Centering and Standardization

## Usage

``` r
scale(x, center = TRUE, scale = TRUE, ...)

# Default S3 method
scale(x, center = TRUE, scale = TRUE, ...)

# S3 method for class 'qts_sample'
scale(
  x,
  center = TRUE,
  scale = TRUE,
  by_row = FALSE,
  keep_summary_stats = FALSE,
  ...
)
```

## Arguments

- x:

  An object coercible into a numeric matrix or an object of class
  [`qts_sample`](https://lmjl-alea.github.io/squat/reference/qts_sample.md)
  representing a sample of observed QTS.

- center:

  A boolean specifying whether to center the sample. If set to `FALSE`,
  the original sample is returned, meaning that no standardization is
  performed regardless of whether argument `scale` was set to `TRUE` or
  not. Defaults to `TRUE`.

- scale:

  A boolean specifying whether to standardize the sample once it has
  been centered. Defaults to `TRUE`.

- ...:

  Extra arguments passed on to next methods.

- by_row:

  A boolean specifying whether the QTS scaling should happen for each
  data point (`by_row = TRUE`) or for each time point
  (`by_row = FALSE`). Defaults to `FALSE`.

- keep_summary_stats:

  A boolean specifying whether the mean and standard deviation used for
  standardizing the data should be stored in the output object. Defaults
  to `FALSE` in which case only the list of properly rescaled QTS is
  returned.

## Value

A list of properly rescaled QTS stored as an object of class
[`qts_sample`](https://lmjl-alea.github.io/squat/reference/qts_sample.md)
when `keep_summary_stats = FALSE`. Otherwise a list with three
components:

- `rescaled_sample`: a list of properly rescaled QTS stored as an object
  of class
  [`qts_sample`](https://lmjl-alea.github.io/squat/reference/qts_sample.md);

- `mean`: a list of numeric vectors storing the corresponding quaternion
  Fréchet means;

- `sd`: a numeric vector storing the corresponding quaternion Fréchet
  standard deviations.

## Examples

``` r
x <- scale(vespa64$igp)
x[[1]]
#> # A tibble: 101 × 5
#>     time         w         x         y         z
#>    <int> <dec:.5!> <dec:.5!> <dec:.5!> <dec:.5!>
#>  1     0   0.99956  -0.01849   0.02258  -0.00497
#>  2     1   0.99960  -0.01996   0.01945  -0.00491
#>  3     2   0.99962  -0.02170   0.01592  -0.00514
#>  4     3   0.99962  -0.02412   0.01211  -0.00573
#>  5     4   0.99960  -0.02614   0.00892  -0.00643
#>  6     5   0.99956  -0.02801   0.00619  -0.00707
#>  7     6   0.99953  -0.02946   0.00407  -0.00758
#>  8     7   0.99951  -0.03028   0.00242  -0.00807
#>  9     8   0.99949  -0.03072   0.00100  -0.00865
#> 10     9   0.99948  -0.03085  -0.00056  -0.00936
#> # ℹ 91 more rows
```
