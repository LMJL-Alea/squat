# PCA for QTS Sample

This is the `S3` specialization of the function
[`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html) for QTS
samples.

## Usage

``` r
# S3 method for class 'qts_sample'
prcomp(x, M = 5, fit = FALSE, ...)
```

## Arguments

- x:

  An object of class
  [qts_sample](https://lmjl-alea.github.io/squat/reference/qts_sample.md).

- M:

  An integer value specifying the number of principal component to
  compute. Defaults to `5L`.

- fit:

  A boolean specifying whether the resulting `prcomp_qts` object should
  store a reconstruction of the sample from the retained PCs. Defaults
  to `FALSE`.

- ...:

  Arguments passed to or from other methods.

## Value

An object of class `prcomp_qts` which is a list with the following
components:

- `x`: An object of class
  [`qts_sample`](https://lmjl-alea.github.io/squat/reference/qts_sample.md)
  as provided by the user, possibly resampled;

- `tpca`: An object of class `MFPCAfit` as produced by the function
  [`MFPCA::MFPCA()`](https://rdrr.io/pkg/MFPCA/man/MFPCA.html),

- `var_props`: A numeric vector storing the percentage of variance
  explained by each PC,

- `total_variance`: A numeric value storing the total variance of the
  sample,

- `mean_qts`: An object of class
  [qts](https://lmjl-alea.github.io/squat/reference/qts.md) containing
  the mean QTS (used for centering the QTS sample before projecting it
  to the tangent space),

- `principal_qts`: A list of
  [qts](https://lmjl-alea.github.io/squat/reference/qts.md)s containing
  the required principal components.

## Details

The `mean_qts` component of the resulting object is the QTS used for
centering. It it part of the `prcomp_qts` object because it is needed to
reconstruct the sample from the retained PCs. The `prcomp_qts` object
also contains the total variance of the sample and the percentage of
variance explained by each PC.

## Examples

``` r
res_pca <- prcomp(vespa64$igp[1:16])
#> ℹ The maximum number of principal component is 15.
```
