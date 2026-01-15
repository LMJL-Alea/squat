# QTS K-Means Alignment Algorithm

This function massages the input quaternion time series to feed them
into the k-means alignment algorithm for jointly clustering and aligning
the input QTS.

## Usage

``` r
kmeans(x, n_clusters, ...)

# Default S3 method
kmeans(
  x,
  n_clusters = 1,
  iter_max = 10,
  nstart = 1,
  algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"),
  trace = FALSE,
  ...
)

# S3 method for class 'qts_sample'
kmeans(
  x,
  n_clusters = 1L,
  seeds = NULL,
  seeding_strategy = c("kmeans++", "exhaustive-kmeans++", "exhaustive", "hclust"),
  is_domain_interval = FALSE,
  transformation = c("identity", "srvf"),
  warping_class = c("none", "shift", "dilation", "affine", "bpd"),
  centroid_type = "mean",
  metric = c("l2", "normalized_l2", "pearson"),
  cluster_on_phase = FALSE,
  use_fence = FALSE,
  ...
)
```

## Arguments

- x:

  Either a numeric matrix of data, or an object that can be coerced to
  such a matrix (such as a numeric vector or a data frame with all
  numeric columns) or an object of class
  [qts_sample](https://lmjl-alea.github.io/squat/reference/qts_sample.md).

- n_clusters:

  An integer value specifying the number of clusters to be look for.

- ...:

  not used.

- iter_max:

  An integer value specifying the maximum number of iterations for
  terminating the k-mean algorithm. Defaults to `10L`.

- nstart:

  if `centers` is a number, how many random sets should be chosen?

- algorithm:

  character: may be abbreviated. Note that `"Lloyd"` and `"Forgy"` are
  alternative names for one algorithm.

- trace:

  logical or integer number, currently only used in the default method
  (`"Hartigan-Wong"`): if positive (or true), tracing information on the
  progress of the algorithm is produced. Higher values may produce more
  tracing information.

- seeds:

  An integer value or vector specifying the indices of the initial
  centroids. If an integer vector, it is interpreted as the indices of
  the intial centroids and should therefore be of length `n_clusters`.
  If an integer value, it is interpreted as the index of the first
  initial centroid and subsequent centroids are chosen according to the
  k-means++ strategy. It can be `NULL` in which case the argument
  `seeding_strategy` is used to automatically provide suitable indices.
  Defaults to `NULL`.

- seeding_strategy:

  A character string specifying the strategy for choosing the initial
  centroids in case the argument `seeds` is set to `NULL`. Choices are
  [`"kmeans++"`](https://en.wikipedia.org/wiki/K-means%2B%2B),
  `"exhaustive-kmeans++"` which performs an exhaustive search over the
  choice of the first centroid, `"exhaustive"` which tries on all
  combinations of initial centroids or `"hclust"` which first performs
  hierarchical clustering using Ward's linkage criterion to identify
  initial centroids. Defaults to `"kmeans++"`, which is the fastest
  strategy.

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

- centroid_type:

  A string specifying the type of centroid to compute. Choices are
  `"mean"`, `"median"` `"medoid"`, `"lowess"` or `"poly"`. Defaults to
  `"mean"`. If LOWESS appproximation is chosen, the user can append an
  integer between 0 and 100 as in `"lowess20"`. This number will be used
  as the smoother span. This gives the proportion of points in the plot
  which influence the smooth at each value. Larger values give more
  smoothness. The default value is 10%. If polynomial approximation is
  chosen, the user can append an positive integer as in `"poly3"`. This
  number will be used as the degree of the polynomial model. The default
  value is `4L`.

- metric:

  A string specifying the metric used to compare curves. Choices are
  `"l2"`, `"normalized_l2"` or `"pearson"`. If
  `transformation == "srvf"`, the metric **must be** `"l2"` because the
  SRVF transform maps absolutely continuous functions to
  square-integrable functions. If `transformation == "identity"` and
  `warping_class` is either `dilation` or `affine`, the metric cab be
  either `"normalized_l2"` or `"pearson"`. The L2 distance is indeed
  **not** dilation-invariant or affine-invariant. The metric can also be
  `"l2"` if `warping_class == "shift"`. Defaults to `"l2"`.

- cluster_on_phase:

  A boolean specifying whether clustering should be based on phase
  variation or amplitude variation. Defaults to `FALSE` which implies
  amplitude variation.

- use_fence:

  A boolean specifying whether the fence algorithm should be used to
  robustify the algorithm against outliers. Defaults to `FALSE`. This is
  used only when `warping_class != "srvf"`.

## Value

An object of class
[`stats::kmeans`](https://rdrr.io/r/stats/kmeans.html) or
[`stats::hclust`](https://rdrr.io/r/stats/hclust.html) or `dbscan_fast`
if the input `x` is NOT of class
[`qts_sample`](https://lmjl-alea.github.io/squat/reference/qts_sample.md).
Otherwise, an object of class `qtsclust` which is effectively a list
with four components:

- `qts_aligned`: An object of class
  [`qts_sample`](https://lmjl-alea.github.io/squat/reference/qts_sample.md)
  storing the sample of aligned QTS;

- `qts_centers`: A list of objects of class
  [`qts`](https://lmjl-alea.github.io/squat/reference/qts.md)
  representing the centers of the clusters;

- `best_clustering`: An object of class
  [`fdacluster::caps`](https://astamm.github.io/fdacluster/reference/caps.html)
  storing the results of the best k-mean alignment result among all
  initialization that were tried.

- `call_name`: A string storing the name of the function that was used
  to produce the clustering structure;

- `call_args`: A list containing the exact arguments that were passed to
  the function `call_name` that produced this output.

## Examples

``` r
out <- kmeans(vespa64$igp[1:10], n_clusters = 2)
```
