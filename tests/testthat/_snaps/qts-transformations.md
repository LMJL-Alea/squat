# the function qts2dts() works

    Code
      qts2dts(vespa64$igp[[1]], vespa64$igp[[2]])
    Output
      # A tibble: 101 x 2
          time distance
         <int>    <dbl>
       1     0  0.0120 
       2     1  0.0106 
       3     2  0.00983
       4     3  0.00936
       5     4  0.00794
       6     5  0.00714
       7     6  0.00731
       8     7  0.00775
       9     8  0.00808
      10     9  0.00824
      # i 91 more rows

# the function qts2nts() works

    Code
      qts2nts(vespa64$igp[[1]], disable_normalization = FALSE)
    Output
      # A tibble: 101 x 2
          time  norm
         <int> <dbl>
       1     0 0.214
       2     1 0.203
       3     2 0.191
       4     3 0.178
       5     4 0.167
       6     5 0.157
       7     6 0.147
       8     7 0.140
       9     8 0.132
      10     9 0.125
      # i 91 more rows

---

    Code
      qts2nts(vespa64$igp[[1]], disable_normalization = TRUE)
    Output
      # A tibble: 101 x 2
          time  norm
         <int> <dbl>
       1     0 0.214
       2     1 0.203
       3     2 0.191
       4     3 0.178
       5     4 0.167
       6     5 0.157
       7     6 0.147
       8     7 0.140
       9     8 0.132
      10     9 0.125
      # i 91 more rows

# the function qts2ats() works

    Code
      qts2ats(vespa64$igp[[1]], disable_normalization = FALSE)
    Output
      # A tibble: 101 x 2
          time  angle
         <int>  <dbl>
       1     0 0     
       2     1 0.0113
       3     2 0.0235
       4     3 0.0366
       5     4 0.0480
       6     5 0.0583
       7     6 0.0673
       8     7 0.0752
       9     8 0.0828
      10     9 0.0908
      # i 91 more rows

---

    Code
      qts2ats(vespa64$igp[[1]], disable_normalization = TRUE)
    Output
      # A tibble: 101 x 2
          time  angle
         <int>  <dbl>
       1     0 0     
       2     1 0.0113
       3     2 0.0235
       4     3 0.0366
       5     4 0.0480
       6     5 0.0583
       7     6 0.0673
       8     7 0.0752
       9     8 0.0828
      10     9 0.0908
      # i 91 more rows

# the function qts2avvts() works

    Code
      qts2avvts(vespa64$igp[[1]])
    Output
      # A tibble: 101 x 4
          time        x        y         z
         <int>    <dbl>    <dbl>     <dbl>
       1     0 -0.0109  -0.00444 -0.000474
       2     1 -0.0103  -0.00519 -0.000562
       3     2 -0.0109  -0.00700 -0.00129 
       4     3 -0.00994 -0.00728 -0.00201 
       5     4 -0.00832 -0.00690 -0.00255 
       6     5 -0.00674 -0.00664 -0.00298 
       7     6 -0.00514 -0.00640 -0.00349 
       8     7 -0.00383 -0.00642 -0.00405 
       9     8 -0.00313 -0.00692 -0.00463 
      10     9 -0.00324 -0.00747 -0.00496 
      # i 91 more rows

---

    Code
      qts2avmts(vespa64$igp[[1]])
    Output
      # A tibble: 101 x 2
          time magnitude
         <int>     <dbl>
       1     0   0.0117 
       2     1   0.0116 
       3     2   0.0130 
       4     3   0.0125 
       5     4   0.0111 
       6     5   0.00992
       7     6   0.00892
       8     7   0.00850
       9     8   0.00889
      10     9   0.00953
      # i 91 more rows

# the function qts2aavts() works

    Code
      qts2aavts(vespa64$igp[[1]])
    Output
      # A tibble: 101 x 4
          time         x         y          z
         <int>     <dbl>     <dbl>      <dbl>
       1     0  0.00234   0.00106   0.000561 
       2     1 -0.00128  -0.00256  -0.000736 
       3     2  0.000199 -0.00104  -0.000725 
       4     3  0.00168   0.000476 -0.000713 
       5     4  0.00156   0.000278 -0.000373 
       6     5  0.00159   0.000252 -0.000467 
       7     6  0.00162   0.000225 -0.000560 
       8     7  0.00100  -0.000258 -0.000572 
       9     8  0.000387 -0.000741 -0.000583 
      10     9 -0.000597 -0.000377 -0.0000683
      # i 91 more rows

---

    Code
      qts2aamts(vespa64$igp[[1]])
    Output
      # A tibble: 101 x 2
          time magnitude
         <int>     <dbl>
       1     0  0.00263 
       2     1  0.00296 
       3     2  0.00129 
       4     3  0.00189 
       5     4  0.00163 
       6     5  0.00168 
       7     6  0.00173 
       8     7  0.00118 
       9     8  0.00102 
      10     9  0.000709
      # i 91 more rows

# the function qts2sqts() works

    Code
      qts2sqts(vespa64$igp[[1]])
    Output
      # A tibble: 101 x 5
          time         w         x         y         z
         <int> <dec:.5!> <dec:.5!> <dec:.5!> <dec:.5!>
       1     0   0.99044   0.07962   0.07962   0.07962
       2     1   0.99163   0.07452   0.07452   0.07452
       3     2   0.99279   0.06920   0.06920   0.06920
       4     3   0.99384   0.06397   0.06397   0.06397
       5     4   0.99469   0.05944   0.05944   0.05944
       6     5   0.99534   0.05569   0.05569   0.05569
       7     6   0.99582   0.05273   0.05273   0.05273
       8     7   0.99616   0.05053   0.05053   0.05053
       9     8   0.99642   0.04882   0.04882   0.04882
      10     9   0.99664   0.04727   0.04727   0.04727
      # i 91 more rows

