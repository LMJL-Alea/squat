# Plot for `qtsclust` objects

This function creates a visualization of the clustering results obtained
on a sample of QTS **without** returning the plot data as an object.

## Usage

``` r
# S3 method for class 'qtsclust'
plot(x, ...)
```

## Arguments

- x:

  An object of class `qtsclust` as produced by
  [`kmeans.qts_sample()`](https://lmjl-alea.github.io/squat/reference/kmeans.md)
  or
  [`hclust.qts_sample()`](https://lmjl-alea.github.io/squat/reference/hclust.md).

- ...:

  Further arguments to be passed to other methods.

## Value

No return value, called for side effects.

## Examples

``` r
out <- kmeans(vespa64$igp[1:10], n_clusters = 2)
plot(out)
```
