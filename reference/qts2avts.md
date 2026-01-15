# QTS Transformation to Angular Velocity Time Series

This function projects a quaternion time series into the space of
angular velocities.

## Usage

``` r
qts2avts(x, body_frame = FALSE)
```

## Arguments

- x:

  An object of class
  [qts](https://lmjl-alea.github.io/squat/reference/qts.md).

- body_frame:

  A boolean specifying whether the fixed frame with respect to which
  coordinates of the angular velocity should be computed is the body
  frame or the global frame. Defaults to `FALSE`.

## Value

A time series stored as a
[tibble::tibble](https://tibble.tidyverse.org/reference/tibble.html)
with columns `time`, `x`, `y` and `z` containing the angular velocity at
each time point.

## Examples

``` r
qts2avts(vespa64$igp[[1]])
#> # A tibble: 100 × 4
#>     time        x        y          z
#>    <dbl>    <dbl>    <dbl>      <dbl>
#>  1     1 -0.0103  -0.00465 -0.0000691
#>  2     2 -0.0105  -0.00624 -0.000639 
#>  3     3 -0.0107  -0.00738 -0.00151  
#>  4     4 -0.00905 -0.00697 -0.00218  
#>  5     5 -0.00766 -0.00673 -0.00274  
#>  6     6 -0.00608 -0.00637 -0.00325  
#>  7     7 -0.00453 -0.00621 -0.00387  
#>  8     8 -0.00365 -0.00642 -0.00451  
#>  9     9 -0.00327 -0.00697 -0.00501  
#> 10    10 -0.00376 -0.00742 -0.00521  
#> # ℹ 90 more rows
```
