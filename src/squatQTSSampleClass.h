#ifndef SQUATQTSSAMPLE_H
#define SQUATQTSSAMPLE_H

#include <Rcpp.h>

//' QTS Geometric Mean
//'
//' This function computes the pointwise geometric mean of a list of QTS.
//'
//' @param qts_list A list of quaternion time series stored as a
//'   \code{\link[tibble]{tibble}}s with columns `time`, `w`, `x`, `y` and `z`.
//'
//' @return A quaternion time series stored as a \code{\link[tibble]{tibble}}
//'   with columns `time`, `w`, `x`, `y` and `z` in which quaternions are the
//'   pointwise geometric mean.
//'
//' @export
//' @examples
//' # TO DO
// [[Rcpp::export]]
Rcpp::DataFrame mean_qts(const Rcpp::List &qts_list);

//' QTS Geometric Median
//'
//' This function computes the pointwise geometric median of a list of QTS.
//'
//' @param qts_list A list of quaternion time series stored as a
//'   \code{\link[tibble]{tibble}}s with columns `time`, `w`, `x`, `y` and `z`.
//'
//' @return A quaternion time series stored as a \code{\link[tibble]{tibble}}
//'   with columns `time`, `w`, `x`, `y` and `z` in which quaternions are the
//'   pointwise geometric median.
//'
//' @export
//' @examples
//' # TO DO
// [[Rcpp::export]]
Rcpp::DataFrame median_qts(const Rcpp::List &qts_list);

#endif /* SQUATQTSSAMPLE_H */
