# squat 0.1.0.9000

# squat 0.1.0

## Major statistical features
* A first API proposal with a class `qts` and a class `qts_sample` for which a
number of methods are properly implemented.
* Available statistical methods for QTS samples: 
  * random generation according to the Gaussian functional model via
  [`rnorm_qts()`](https://lmjl-alea.github.io/squat/reference/rnorm_qts.html),
  * [`scale()`](https://lmjl-alea.github.io/squat/reference/scale.html), 
  * [`mean()`](https://lmjl-alea.github.io/squat/reference/mean.qts_sample.html), 
  * [`median()`](https://lmjl-alea.github.io/squat/reference/median.qts_sample.html),
  * distance matrix computation via [`distDTW()`](https://lmjl-alea.github.io/squat/reference/distDTW.html) (i.e. for now we use the dynamic time warping),
  * tangent principal component analysis via [`prcomp()`](https://lmjl-alea.github.io/squat/reference/prcomp.qts_sample.html),
  * k-means with optional alignment via [`kmeans()`](https://lmjl-alea.github.io/squat/reference/kmeans.html).
* Added multiple ways of displaying samples of QTS.
* Added two example datasets.

## Improvements
* Make all functions applicable to a single QTS also applicable to QTS samples,
with appropriate class for the output.
* Enable `as_qts_sample()` to generate a QTS sample of size 1 from a single QTS
as input argument.
* Rename `change_points` argument to the `plot.qts()` function to better reflect
its flexibility.
* Added subset operator for QTS sample objects.
* Added `append` S3 method for QTS sample objects.
* Added `hemispherize()` function to remove any discontinuities in QTS due to
quaternion flips.
* Any parallelization computation is now handled using the
[**futureverse**](https://www.futureverse.org) principles and, in particular,
implemented through the use of the [**furrr**](https://furrr.futureverse.org)
package.

# squat 0.0.1

* Added a `NEWS.md` file to track changes to the package.
* Initial version.
