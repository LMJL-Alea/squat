# The function derivative_qts() works

    Code
      derivative_qts(vespa64$igp[[1]])
    Output
      # A tibble: 100 x 5
          time     w        x        y         z
         <int> <dbl>    <dbl>    <dbl>     <dbl>
       1     1  1.00 -0.00519 -0.00222 -0.000398
       2     2  1.00 -0.00529 -0.00305 -0.000574
       3     3  1.00 -0.00534 -0.00367 -0.000939
       4     4  1.00 -0.00447 -0.00352 -0.00119 
       5     5  1.00 -0.00374 -0.00345 -0.00140 
       6     6  1.00 -0.00292 -0.00331 -0.00158 
       7     7  1.00 -0.00212 -0.00327 -0.00182 
       8     8  1.00 -0.00166 -0.00341 -0.00209 
       9     9  1.00 -0.00145 -0.00371 -0.00229 
      10    10  1.00 -0.00168 -0.00394 -0.00239 
      # ... with 90 more rows

# Logarithm and exponential for QTS work

    Code
      x
    Output
      # A tibble: 101 x 5
          time         w      x      y       z
         <int>     <dbl>  <dbl>  <dbl>   <dbl>
       1     0 -1.11e-16 0.0799 0.0700 0.0134 
       2     1  0        0.0747 0.0677 0.0132 
       3     2  2.22e-16 0.0694 0.0647 0.0127 
       4     3  0        0.0641 0.0610 0.0119 
       5     4  0        0.0596 0.0575 0.0107 
       6     5  0        0.0558 0.0541 0.00933
       7     6  0        0.0528 0.0508 0.00772
       8     7 -1.11e-16 0.0506 0.0476 0.00584
       9     8 -1.11e-16 0.0489 0.0443 0.00366
      10     9  0        0.0473 0.0407 0.00126
      # ... with 91 more rows

# Function reorient_qts() works

    Code
      reorient_qts(vespa64$igp[[1]], disable_normalization = FALSE)
    Output
      # A tibble: 101 x 5
          time     w         x         y         z
         <int> <dbl>     <dbl>     <dbl>     <dbl>
       1     0 1     -1.39e-17  2.17e-19 -8.67e-19
       2     1 1.00  -5.19e- 3 -2.22e- 3 -3.98e- 4
       3     2 1.00  -1.05e- 2 -5.27e- 3 -9.68e- 4
       4     3 1.00  -1.58e- 2 -8.95e- 3 -1.90e- 3
       5     4 1.00  -2.03e- 2 -1.25e- 2 -3.07e- 3
       6     5 1.00  -2.40e- 2 -1.59e- 2 -4.44e- 3
       7     6 0.999 -2.69e- 2 -1.93e- 2 -5.99e- 3
       8     7 0.999 -2.90e- 2 -2.26e- 2 -7.77e- 3
       9     8 0.999 -3.07e- 2 -2.60e- 2 -9.79e- 3
      10     9 0.999 -3.21e- 2 -2.98e- 2 -1.20e- 2
      # ... with 91 more rows

---

    Code
      reorient_qts(vespa64$igp[[1]], disable_normalization = TRUE)
    Output
      # A tibble: 101 x 5
          time     w         x        y         z
         <int> <dbl>     <dbl>    <dbl>     <dbl>
       1     0 1     -1.08e-19  0        8.67e-19
       2     1 1.00  -5.19e- 3 -0.00222 -3.98e- 4
       3     2 1.00  -1.05e- 2 -0.00527 -9.68e- 4
       4     3 1.00  -1.58e- 2 -0.00895 -1.90e- 3
       5     4 1.00  -2.03e- 2 -0.0125  -3.07e- 3
       6     5 1.00  -2.40e- 2 -0.0159  -4.44e- 3
       7     6 0.999 -2.69e- 2 -0.0193  -5.99e- 3
       8     7 0.999 -2.90e- 2 -0.0226  -7.77e- 3
       9     8 0.999 -3.07e- 2 -0.0260  -9.79e- 3
      10     9 0.999 -3.21e- 2 -0.0298  -1.20e- 2
      # ... with 91 more rows

# Function normalize_qts() works

    Code
      normalize_qts(vespa64$igp[[1]])
    Output
      # A tibble: 101 x 5
          time     w      x      y       z
         <int> <dbl>  <dbl>  <dbl>   <dbl>
       1     0 0.994 0.0797 0.0699 0.0133 
       2     1 0.995 0.0746 0.0676 0.0131 
       3     2 0.995 0.0693 0.0646 0.0127 
       4     3 0.996 0.0640 0.0609 0.0118 
       5     4 0.997 0.0595 0.0574 0.0107 
       6     5 0.997 0.0557 0.0540 0.00932
       7     6 0.997 0.0527 0.0508 0.00772
       8     7 0.998 0.0506 0.0476 0.00583
       9     8 0.998 0.0488 0.0443 0.00366
      10     9 0.998 0.0473 0.0407 0.00125
      # ... with 91 more rows

# Function centring_qts() works (standardize = FALSE, keep_summary_stats = FALSE)

    Code
      centring_qts(x = vespa64$igp[[1]], standardize = FALSE, keep_summary_stats = FALSE)
    Output
      # A tibble: 101 x 5
          time     w      x      y       z
         <int> <dbl>  <dbl>  <dbl>   <dbl>
       1     0 0.994 0.0938 0.0605 0.0148 
       2     1 0.994 0.0887 0.0582 0.0145 
       3     2 0.995 0.0834 0.0552 0.0140 
       4     3 0.996 0.0781 0.0515 0.0131 
       5     4 0.996 0.0736 0.0481 0.0118 
       6     5 0.996 0.0699 0.0447 0.0104 
       7     6 0.997 0.0669 0.0415 0.00869
       8     7 0.997 0.0648 0.0383 0.00674
       9     8 0.997 0.0631 0.0350 0.00451
      10     9 0.998 0.0616 0.0315 0.00204
      # ... with 91 more rows

# Function centring_qts() works (standardize = TRUE, keep_summary_stats = FALSE)

    Code
      centring_qts(x = vespa64$igp[[1]], standardize = TRUE, keep_summary_stats = FALSE)
    Output
      # A tibble: 101 x 5
          time     w      x      y       z
         <int> <dbl>  <dbl>  <dbl>   <dbl>
       1     0 0.997 0.0624 0.0402 0.00987
       2     1 0.997 0.0590 0.0387 0.00967
       3     2 0.998 0.0555 0.0367 0.00932
       4     3 0.998 0.0519 0.0343 0.00869
       5     4 0.998 0.0490 0.0319 0.00787
       6     5 0.998 0.0465 0.0297 0.00689
       7     6 0.999 0.0445 0.0276 0.00578
       8     7 0.999 0.0431 0.0255 0.00448
       9     8 0.999 0.0419 0.0233 0.00300
      10     9 0.999 0.0409 0.0209 0.00135
      # ... with 91 more rows

# Function centring_qts() works (standardize = FALSE, keep_summary_stats = TRUE)

    Code
      centring_qts(x = vespa64$igp[[1]], standardize = FALSE, keep_summary_stats = TRUE)
    Output
      $qts
      # A tibble: 101 x 5
          time     w      x      y       z
         <int> <dbl>  <dbl>  <dbl>   <dbl>
       1     0 0.994 0.0938 0.0605 0.0148 
       2     1 0.994 0.0887 0.0582 0.0145 
       3     2 0.995 0.0834 0.0552 0.0140 
       4     3 0.996 0.0781 0.0515 0.0131 
       5     4 0.996 0.0736 0.0481 0.0118 
       6     5 0.996 0.0699 0.0447 0.0104 
       7     6 0.997 0.0669 0.0415 0.00869
       8     7 0.997 0.0648 0.0383 0.00674
       9     8 0.997 0.0631 0.0350 0.00451
      10     9 0.998 0.0616 0.0315 0.00204
      # ... with 91 more rows
      
      $mean
      [1]  0.9998551535 -0.0143000963  0.0092262363  0.0002360886
      
      $sd
      [1] 0
      

