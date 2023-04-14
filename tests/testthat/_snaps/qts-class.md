# Function centring() works (standardize = FALSE, keep_summary_stats = FALSE)

    Code
      centring(x = vespa64$igp[[1]], standardize = FALSE, keep_summary_stats = FALSE)
    Output
      # A tibble: 101 x 5
          time         w         x         y         z
         <int> <dec:.5!> <dec:.5!> <dec:.5!> <dec:.5!>
       1     0   0.99364   0.09383   0.06049   0.01484
       2     1   0.99425   0.08868   0.05824   0.01455
       3     2   0.99489   0.08343   0.05518   0.01401
       4     3   0.99553   0.07811   0.05153   0.01307
       5     4   0.99606   0.07365   0.04805   0.01184
       6     5   0.99650   0.06989   0.04468   0.01037
       7     6   0.99686   0.06694   0.04146   0.00869
       8     7   0.99714   0.06478   0.03830   0.00674
       9     8   0.99738   0.06307   0.03503   0.00451
      10     9   0.99761   0.06156   0.03147   0.00204
      # ... with 91 more rows

# Function centring() works (standardize = TRUE, keep_summary_stats = FALSE)

    Code
      centring(x = vespa64$igp[[1]], standardize = TRUE, keep_summary_stats = FALSE)
    Output
      # A tibble: 101 x 5
          time         w         x         y         z
         <int> <dec:.5!> <dec:.5!> <dec:.5!> <dec:.5!>
       1     0   0.99719   0.06241   0.04023   0.00987
       2     1   0.99746   0.05897   0.03873   0.00967
       3     2   0.99774   0.05548   0.03669   0.00932
       4     3   0.99802   0.05194   0.03426   0.00869
       5     4   0.99826   0.04896   0.03195   0.00787
       6     5   0.99845   0.04646   0.02970   0.00689
       7     6   0.99861   0.04450   0.02756   0.00578
       8     7   0.99874   0.04306   0.02546   0.00448
       9     8   0.99885   0.04192   0.02328   0.00300
      10     9   0.99894   0.04091   0.02091   0.00135
      # ... with 91 more rows

# Function centring() works (standardize = FALSE, keep_summary_stats = TRUE)

    Code
      centring(x = vespa64$igp[[1]], standardize = FALSE, keep_summary_stats = TRUE)
    Output
      $qts
      # A tibble: 101 x 5
          time         w         x         y         z
         <int> <dec:.5!> <dec:.5!> <dec:.5!> <dec:.5!>
       1     0   0.99364   0.09383   0.06049   0.01484
       2     1   0.99425   0.08868   0.05824   0.01455
       3     2   0.99489   0.08343   0.05518   0.01401
       4     3   0.99553   0.07811   0.05153   0.01307
       5     4   0.99606   0.07365   0.04805   0.01184
       6     5   0.99650   0.06989   0.04468   0.01037
       7     6   0.99686   0.06694   0.04146   0.00869
       8     7   0.99714   0.06478   0.03830   0.00674
       9     8   0.99738   0.06307   0.03503   0.00451
      10     9   0.99761   0.06156   0.03147   0.00204
      # ... with 91 more rows
      
      $mean
      [1]  0.9998551535 -0.0143000963  0.0092262363  0.0002360886
      
      $sd
      [1] 0
      
