# Dynamic Time Warping for Quaternion Time Series

This function evaluates the Dynamic Time Warping (DTW) distance between
two quaternion time series (QTS).

## Usage

``` r
DTW(
  qts1,
  qts2,
  resample = TRUE,
  disable_normalization = FALSE,
  distance_only = FALSE,
  step_pattern = dtw::symmetric2
)
```

## Arguments

- qts1:

  An object of class
  [qts](https://lmjl-alea.github.io/squat/reference/qts.md).

- qts2:

  An object of class
  [qts](https://lmjl-alea.github.io/squat/reference/qts.md).

- resample:

  A boolean specifying whether the QTS should be uniformly resampled on
  their domain before computing distances. Defaults to `TRUE`.

- disable_normalization:

  A boolean specifying whether quaternion normalization should be
  disabled. Defaults to `FALSE` which ensures that we always deal with
  unit quaternions.

- distance_only:

  A boolean specifying whether to only compute distance (no backtrack,
  faster). Defaults to `FALSE`.

- step_pattern:

  A [dtw::stepPattern](https://rdrr.io/pkg/dtw/man/stepPattern.html)
  specifying the local constraints on the warping path. Defaults to
  [dtw::symmetric2](https://rdrr.io/pkg/dtw/man/stepPattern.html) which
  uses symmetric and normalizable warping paths with no local slope
  constraints. See
  [dtw::stepPattern](https://rdrr.io/pkg/dtw/man/stepPattern.html) for
  more information.

## Value

An object of class [dtw::dtw](https://rdrr.io/pkg/dtw/man/dtw.html)
storing the dynamic time warping results.

## Details

If no evaluation grid is provided, the function assumes that the two
input QTS are evaluated on the same grid.

## Examples

``` r
DTW(vespa64$igp[[1]], vespa64$igp[[2]])
#> DTW alignment object
#>  Alignment size (query x reference): 101 x 101
#>  Call: dtw::dtw(x = M, step.pattern = step_pattern, distance.only = distance_only)
```
