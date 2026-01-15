# Plot for `prcomp_qts` objects

This function creates a visualization of the results of the PCA applied
on a sample of QTS **without** returning the plot data as an object.

## Usage

``` r
# S3 method for class 'prcomp_qts'
plot(x, what = "PC1", ...)

# S3 method for class 'prcomp_qts'
screeplot(x, ...)
```

## Arguments

- x:

  An object of class `prcomp_qts` as produced by the
  [`prcomp.qts_sample()`](https://lmjl-alea.github.io/squat/reference/prcomp.qts_sample.md)
  method.

- what:

  A string specifying what kind of visualization the user wants to
  perform. Choices are words starting with `PC` and ending with a PC
  number (in which case the mean QTS is displayed along with its
  perturbations due to the required PC) or `scores` (in which case
  individuals are projected on the required plane). Defaults to `PC1`.

- ...:

  If `what = "PC?"`, the user can specify whether to plot the QTS in the
  tangent space or in the original space by providing a boolean argument
  `original_space` which defaults to `TRUE`. If `what = "scores"`, the
  user can specify the plane onto which the individuals will be
  projected by providing a length-2 integer vector argument `plane`
  which defaults to `1:2`.

## Value

No return value, called for side effects.

## Examples

``` r
df <- as_qts_sample(vespa64$igp[1:16])
res_pca <- prcomp(df)
#> ℹ The maximum number of principal component is 15.

# You can plot the effect of a PC on the mean
plot(res_pca, what = "PC1")
#> The `original_space` boolean argument is not specified. Defaulting to TRUE.


# You can plot the data points in a PC plane
plot(res_pca, what = "scores")
#> The `plane` length-2 integer vector argument is not specified. Defaulting to
#> 1:2.
```
