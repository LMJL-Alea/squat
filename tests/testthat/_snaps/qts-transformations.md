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
      # ... with 91 more rows

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
      # ... with 91 more rows

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
      # ... with 91 more rows

