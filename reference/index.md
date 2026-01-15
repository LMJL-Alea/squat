# Package index

## QTS Class

- [`as_qts()`](https://lmjl-alea.github.io/squat/reference/qts.md)
  [`is_qts()`](https://lmjl-alea.github.io/squat/reference/qts.md)
  [`format(`*`<qts>`*`)`](https://lmjl-alea.github.io/squat/reference/qts.md)
  [`print(`*`<qts>`*`)`](https://lmjl-alea.github.io/squat/reference/qts.md)
  : QTS Class

- [`centring()`](https://lmjl-alea.github.io/squat/reference/centring.md)
  : QTS Centering and Standardization

- [`autoplot(`*`<qts>`*`)`](https://lmjl-alea.github.io/squat/reference/autoplot.qts.md)
  :

  Plot for `qts` objects

- [`plot(`*`<qts>`*`)`](https://lmjl-alea.github.io/squat/reference/plot.qts.md)
  :

  Plot for `qts` objects

- [`` `+`( ``*`<qts>`*`)`](https://lmjl-alea.github.io/squat/reference/plus-.qts.md)
  :

  Operator `+` for `qts` Objects

- [`` `-`( ``*`<qts>`*`)`](https://lmjl-alea.github.io/squat/reference/dot-qts.md)
  :

  Operator `-` for `qts` Objects

- [`` `*`( ``*`<qts>`*`)`](https://lmjl-alea.github.io/squat/reference/times-.qts.md)
  :

  Operator `*` for `qts` Objects

- [`inverse_qts()`](https://lmjl-alea.github.io/squat/reference/inverse_qts.md)
  :

  Inverse Operator for `qts` Objects

## QTS Sample Class

- [`as_qts_sample()`](https://lmjl-alea.github.io/squat/reference/qts_sample.md)
  [`is_qts_sample()`](https://lmjl-alea.github.io/squat/reference/qts_sample.md)
  [`` `[`( ``*`<qts_sample>`*`)`](https://lmjl-alea.github.io/squat/reference/qts_sample.md)
  : QTS Sample Class

- [`append()`](https://lmjl-alea.github.io/squat/reference/append.md) :
  QTS Sample Concatenation

- [`rnorm_qts()`](https://lmjl-alea.github.io/squat/reference/rnorm_qts.md)
  : QTS Random Sampling

- [`scale()`](https://lmjl-alea.github.io/squat/reference/scale.md) :
  QTS Sample Centering and Standardization

- [`mean(`*`<qts_sample>`*`)`](https://lmjl-alea.github.io/squat/reference/mean.qts_sample.md)
  : QTS Geometric Mean

- [`median(`*`<qts_sample>`*`)`](https://lmjl-alea.github.io/squat/reference/median.qts_sample.md)
  : QTS Geometric Median

- [`autoplot(`*`<qts_sample>`*`)`](https://lmjl-alea.github.io/squat/reference/autoplot.qts_sample.md)
  :

  Plot for `qts_sample` objects

- [`plot(`*`<qts_sample>`*`)`](https://lmjl-alea.github.io/squat/reference/plot.qts_sample.md)
  :

  Plot for `qts_sample` objects

## QTS Wrangling

- [`differentiate()`](https://lmjl-alea.github.io/squat/reference/differentiate.md)
  : QTS Differentiation
- [`straighten()`](https://lmjl-alea.github.io/squat/reference/straighten.md)
  : QTS Straightening
- [`log(`*`<qts>`*`)`](https://lmjl-alea.github.io/squat/reference/log.md)
  [`log(`*`<qts_sample>`*`)`](https://lmjl-alea.github.io/squat/reference/log.md)
  : QTS Logarithm
- [`exp(`*`<qts>`*`)`](https://lmjl-alea.github.io/squat/reference/exp.md)
  [`exp(`*`<qts_sample>`*`)`](https://lmjl-alea.github.io/squat/reference/exp.md)
  : QTS Exponential
- [`reorient()`](https://lmjl-alea.github.io/squat/reference/reorient.md)
  : QTS Reorientation
- [`normalize()`](https://lmjl-alea.github.io/squat/reference/normalize.md)
  : QTS Normalization
- [`resample()`](https://lmjl-alea.github.io/squat/reference/resample.md)
  : QTS Resampling
- [`smooth()`](https://lmjl-alea.github.io/squat/reference/smooth.md) :
  QTS Smoothing via SLERP Interpolation
- [`hemispherize()`](https://lmjl-alea.github.io/squat/reference/hemispherize.md)
  : QTS Hemispherization
- [`moving_average()`](https://lmjl-alea.github.io/squat/reference/moving_average.md)
  : QTS Moving Average

## QTS Transformations

- [`qts2dts()`](https://lmjl-alea.github.io/squat/reference/qts2dts.md)
  : QTS Transformation To Distance Time Series
- [`qts2aats()`](https://lmjl-alea.github.io/squat/reference/qts2aats.md)
  : QTS Transformation to Angle-Axis Time Series
- [`qts2nts()`](https://lmjl-alea.github.io/squat/reference/qts2nts.md)
  : QTS Transformation To Norm Time Series
- [`qts2ats()`](https://lmjl-alea.github.io/squat/reference/qts2ats.md)
  : QTS Transformation To Angle Time Series
- [`qts2avts()`](https://lmjl-alea.github.io/squat/reference/qts2avts.md)
  : QTS Transformation to Angular Velocity Time Series
- [`qts2rpyts()`](https://lmjl-alea.github.io/squat/reference/qts2rpyts.md)
  : QTS Transformation to Roll-Pitch-Yaw Time Series

## Principal Component Analysis

- [`prcomp(`*`<qts_sample>`*`)`](https://lmjl-alea.github.io/squat/reference/prcomp.qts_sample.md)
  : PCA for QTS Sample

- [`predict(`*`<prcomp_qts>`*`)`](https://lmjl-alea.github.io/squat/reference/predict.prcomp_qts.md)
  : Predict QTS from PCA decomposition

- [`autoplot(`*`<prcomp_qts>`*`)`](https://lmjl-alea.github.io/squat/reference/autoplot.prcomp_qts.md)
  :

  Plot for `prcomp_qts` objects

- [`plot(`*`<prcomp_qts>`*`)`](https://lmjl-alea.github.io/squat/reference/plot.prcomp_qts.md)
  [`screeplot(`*`<prcomp_qts>`*`)`](https://lmjl-alea.github.io/squat/reference/plot.prcomp_qts.md)
  :

  Plot for `prcomp_qts` objects

## Clustering

- [`kmeans()`](https://lmjl-alea.github.io/squat/reference/kmeans.md) :
  QTS K-Means Alignment Algorithm

- [`hclust()`](https://lmjl-alea.github.io/squat/reference/hclust.md) :
  QTS Hierarchical Agglomerative Clustering

- [`dbscan()`](https://lmjl-alea.github.io/squat/reference/dbscan.md) :
  QTS Nearest-Neighbor Clustering

- [`autoplot(`*`<qtsclust>`*`)`](https://lmjl-alea.github.io/squat/reference/autoplot.qtsclust.md)
  :

  Plot for `qtsclust` objects

- [`plot(`*`<qtsclust>`*`)`](https://lmjl-alea.github.io/squat/reference/plot.qtsclust.md)
  :

  Plot for `qtsclust` objects

## Distances

- [`DTW()`](https://lmjl-alea.github.io/squat/reference/DTW.md) :
  Dynamic Time Warping for Quaternion Time Series
- [`dist()`](https://lmjl-alea.github.io/squat/reference/dist.md) : QTS
  Distance Matrix Computation

## Datasets

- [`vespa`](https://lmjl-alea.github.io/squat/reference/vespa.md) : The
  VESPA dataset
- [`vespa64`](https://lmjl-alea.github.io/squat/reference/vespa64.md) :
  The VESPA64 dataset
