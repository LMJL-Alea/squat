# QTS Sample Class

A collection of functions that implements the QTS sample class. It
currently provides the `as_qts_sample()` function for QTS sample
coercion of lists of
[`qts`](https://lmjl-alea.github.io/squat/reference/qts.md) objects, the
`is_qts_sample()` function for checking if an object is a QTS sample and
the subset operator.

## Usage

``` r
as_qts_sample(x)

is_qts_sample(x)

# S3 method for class 'qts_sample'
x[i, simplify = FALSE]
```

## Arguments

- x:

  A list of
  [`tibble::tibble`](https://tibble.tidyverse.org/reference/tibble.html)s,
  each of which with columns `time`, `w`, `x`, `y` and `z`.

- i:

  A valid expression to subset observations from a QTS sample.

- simplify:

  A boolean value specifying whether the resulting subset should be
  turned into a single QTS in case the subset is of size 1. Defaults to
  `FALSE`.

## Value

An object of class `qts_sample`.

## Details

A QTS sample is a collection of quaternion time series (QTS), each of
which is stored as a
[`tibble::tibble`](https://tibble.tidyverse.org/reference/tibble.html)
with 5 columns:

- `time`: A first column specifying the time points at which quaternions
  were collected;

- `w`: A second column specifying the first coordinate of the collected
  quaternions;

- `x`: A third column specifying the second coordinate of the collected
  quaternions;

- `y`: A fourth column specifying the third coordinate of the collected
  quaternions;

- `z`: A fifth column specifying the fourth coordinate of the collected
  quaternions.

## Examples

``` r
x <- vespa64$igp
y <- as_qts_sample(x)
is_qts_sample(x)
#> [1] TRUE
is_qts_sample(y)
#> [1] TRUE
x[1]
#> [[1]]
#> # A tibble: 101 × 5
#>     time         w         x         y         z
#>    <int> <dec:.5!> <dec:.5!> <dec:.5!> <dec:.5!>
#>  1     0   0.99427   0.07973   0.06988   0.01334
#>  2     1   0.99483   0.07457   0.06763   0.01313
#>  3     2   0.99542   0.06931   0.06457   0.01269
#>  4     3   0.99602   0.06398   0.06091   0.01184
#>  5     4   0.99652   0.05949   0.05742   0.01070
#>  6     5   0.99694   0.05572   0.05403   0.00932
#>  7     6   0.99729   0.05275   0.05079   0.00772
#>  8     7   0.99757   0.05056   0.04761   0.00583
#>  9     8   0.99782   0.04883   0.04430   0.00366
#> 10     9   0.99805   0.04730   0.04071   0.00125
#> # ℹ 91 more rows
#> 
#> attr(,"class")
#> [1] "qts_sample" "list"      
x[1, simplify = TRUE]
#> # A tibble: 101 × 5
#>     time         w         x         y         z
#>    <int> <dec:.5!> <dec:.5!> <dec:.5!> <dec:.5!>
#>  1     0   0.99427   0.07973   0.06988   0.01334
#>  2     1   0.99483   0.07457   0.06763   0.01313
#>  3     2   0.99542   0.06931   0.06457   0.01269
#>  4     3   0.99602   0.06398   0.06091   0.01184
#>  5     4   0.99652   0.05949   0.05742   0.01070
#>  6     5   0.99694   0.05572   0.05403   0.00932
#>  7     6   0.99729   0.05275   0.05079   0.00772
#>  8     7   0.99757   0.05056   0.04761   0.00583
#>  9     8   0.99782   0.04883   0.04430   0.00366
#> 10     9   0.99805   0.04730   0.04071   0.00125
#> # ℹ 91 more rows
```
