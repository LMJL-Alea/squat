# Predict QTS from PCA decomposition

This function predicts the QTS of a new sample from the PCA
decomposition of a previous sample.

## Usage

``` r
# S3 method for class 'prcomp_qts'
predict(object, newdata, ...)
```

## Arguments

- object:

  An object of class `prcomp_qts` as produced by the
  [`prcomp.qts_sample()`](https://lmjl-alea.github.io/squat/reference/prcomp.qts_sample.md)
  method.

- newdata:

  An object of class
  [`qts`](https://lmjl-alea.github.io/squat/reference/qts.md) or
  [`qts_sample`](https://lmjl-alea.github.io/squat/reference/qts_sample.md)
  specifying a QTS or a sample of QTS. The QTS should be evaluated on
  the same grid as the one used to fit the PCA model. If the evaluation
  grids map the same domain but with different sampling frequenciesa,
  the QTS will be linearly interpolated (in the Lie algebra) to the
  common grid used to fit the PCA model.

- ...:

  Additional arguments. Not used here.

## Value

An object of class
[`qts_sample`](https://lmjl-alea.github.io/squat/reference/qts_sample.md)
containing the predicted QTS.

## Examples

``` r
# Fit PCA model
pr <- prcomp(vespa64$igp, M = 5)
#> ℹ The maximum number of principal component is 63.

# Predict QTS
new_qts <- predict(pr)
```
