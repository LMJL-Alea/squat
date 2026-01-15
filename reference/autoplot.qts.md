# Plot for [`qts`](https://lmjl-alea.github.io/squat/reference/qts.md) objects

This function creates a visualization of a QTS and returns the
corresponding
[ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object which enable further customization of the plot.

## Usage

``` r
# S3 method for class 'qts'
autoplot(object, highlighted_points = NULL, ...)
```

## Arguments

- object:

  An object of class
  [qts](https://lmjl-alea.github.io/squat/reference/qts.md).

- highlighted_points:

  An integer vector specifying point indices to be highlighted. Defaults
  to `NULL`, in which case no point will be highlighted with respect to
  the others.

- ...:

  Further arguments to be passed on to next methods.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Examples

``` r
ggplot2::autoplot(vespa64$igp[[1]])
```
