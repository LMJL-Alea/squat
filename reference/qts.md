# QTS Class

A collection of functions that implements the QTS class. It currently
provides the `as_qts()` function for QTS coercion of
[`tibble::tibble`](https://tibble.tidyverse.org/reference/tibble.html)s
and the `is_qts()` function for checking if an object is a QTS.

## Usage

``` r
as_qts(x)

is_qts(x)

# S3 method for class 'qts'
format(x, digits = 5, ...)

# S3 method for class 'qts'
print(x, ...)
```

## Arguments

- x:

  A
  [`tibble::tibble`](https://tibble.tidyverse.org/reference/tibble.html)
  with columns `time`, `w`, `x`, `y` and `z`.

- digits:

  An integer value specifying the number of digits to keep for printing.
  Defaults to `5L`.

- ...:

  Further arguments passed to or from other methods.

## Value

An object of class qts.

## Details

A quaternion time series (QTS) is stored as a
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
qts1 <- vespa64$igp[[1]]
qts2 <- as_qts(qts1)
is_qts(qts1)
#> [1] TRUE
is_qts(qts2)
#> [1] TRUE
```
