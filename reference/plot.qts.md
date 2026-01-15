# Plot for [`qts`](https://lmjl-alea.github.io/squat/reference/qts.md) objects

This function creates a visualization of a QTS **without** returning the
plot data as an object.

## Usage

``` r
# S3 method for class 'qts'
plot(x, highlighted_points = NULL, ...)
```

## Arguments

- x:

  An object of class
  [qts](https://lmjl-alea.github.io/squat/reference/qts.md).

- highlighted_points:

  An integer vector specifying point indices to be highlighted. Defaults
  to `NULL`, in which case no point will be highlighted with respect to
  the others.

- ...:

  Further arguments to be passed on to next methods.

## Value

No return value, called for side effects.

## Examples

``` r
plot(vespa64$igp[[1]])
```
