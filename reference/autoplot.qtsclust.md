# Plot for `qtsclust` objects

This function creates a visualization of the clustering results obtained
on a sample of QTS and returns the corresponding
[ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object which enable further customization of the plot.

## Usage

``` r
# S3 method for class 'qtsclust'
autoplot(object, ...)
```

## Arguments

- object:

  An object of class `qtsclust` as produced by
  [`kmeans.qts_sample()`](https://lmjl-alea.github.io/squat/reference/kmeans.md)
  or
  [`hclust.qts_sample()`](https://lmjl-alea.github.io/squat/reference/hclust.md).

- ...:

  Further arguments to be passed to other methods.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Examples

``` r
out <- kmeans(vespa64$igp[1:10], n_clusters = 2)
ggplot2::autoplot(out)
```
