# QTS Distance Matrix Computation

This function massages an input sample of quaternion time series to turn
it into a pairwise distance matrix.

## Usage

``` r
dist(x, metric, ...)

# Default S3 method
dist(
  x,
  metric = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
  diag = FALSE,
  upper = FALSE,
  p = 2,
  ...
)

# S3 method for class 'qts_sample'
dist(
  x,
  metric = c("l2", "normalized_l2", "pearson", "dtw"),
  is_domain_interval = FALSE,
  transformation = c("identity", "srvf"),
  warping_class = c("none", "shift", "dilation", "affine", "bpd"),
  rotation_invariance = FALSE,
  cluster_on_phase = FALSE,
  labels = NULL,
  ncores = 1L,
  ...
)
```

## Arguments

- x:

  A numeric matrix, data frame,
  [stats::dist](https://rdrr.io/r/stats/dist.html) object or object of
  class
  [qts_sample](https://lmjl-alea.github.io/squat/reference/qts_sample.md)
  specifying the sample on which to compute the pairwise distance
  matrix.

- metric:

  A character string specifying the distance measure to be used. This
  must be one of `"euclidean"`, `"maximum"`, `"manhattan"`,
  `"canberra"`, `"binary"` or `"minkowski"` if `x` is not a QTS sample.
  Otherwise, it must be one of `"l2"`, `"pearson"` or `"dtw"`.

- ...:

  not used.

- diag:

  logical value indicating whether the diagonal of the distance matrix
  should be printed by `print.dist`.

- upper:

  logical value indicating whether the upper triangle of the distance
  matrix should be printed by `print.dist`.

- p:

  The power of the Minkowski distance.

- is_domain_interval:

  A boolean specifying whether the sample of curves is defined on a
  fixed interval. Defaults to `FALSE`.

- transformation:

  A string specifying the transformation to apply to the original sample
  of curves. Choices are no transformation
  (`transformation = "identity"`) or square-root velocity function
  `transformation = "srvf"`. Defaults to `"identity"`.

- warping_class:

  A string specifying the class of warping functions. Choices are no
  warping (`warping_class = "none"`), shift `y = x + b`
  (`warping_class = "shift"`), dilation `y = ax`
  (`warping_class = "dilation"`), affine `y = ax + b`
  (`warping_class = "affine"`) or boundary-preserving diffeomorphism
  (`warping_class = "bpd"`). Defaults to `"none"`.

- rotation_invariance:

  A boolean value specifying whether the distance should be invariant to
  rotation. This is only relevant when `is_domain_interval` is `TRUE`
  and `transformation` is `"srvf"` and `warped_class` is `"bpd"`.
  Defaults to `FALSE`.

- cluster_on_phase:

  A boolean specifying whether clustering should be based on phase
  variation or amplitude variation. Defaults to `FALSE` which implies
  amplitude variation.

- labels:

  A character vector specifying curve labels. Defaults to `NULL` which
  uses sequential numbers as labels.

- ncores:

  An integer value specifying the number of cores to use for parallel
  computation. Defaults to `1`.

## Value

An object of class [stats::dist](https://rdrr.io/r/stats/dist.html).

## Examples

``` r
D <- dist(vespa64$igp[1:5])
```
