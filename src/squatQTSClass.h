#ifndef SQUATQTSCLASS_H
#define SQUATQTSCLASS_H

#include <Rcpp.h>

//' QTS Reorientation
//'
//' This function reorients the quaternions in a QTS for representing attitude
//' with respect to the first time point.
//'
//' @param qts A quaternion time series stored as a \code{\link[tibble]{tibble}}
//'   with columns `time`, `w`, `x`, `y` and `z`.
//' @param disable_normalization A boolean specifying whether quaternion
//'   normalization should be disabled. Defaults to `FALSE`.
//'
//' @return A quaternion time series stored as a \code{\link[tibble]{tibble}}
//'   with columns `time`, `w`, `x`, `y` and `z` in which quaternions measure
//'   attitude with respect to the first time point.
//'
//' @export
//' @examples
//' # TO DO
// [[Rcpp::export]]
Rcpp::DataFrame reorient_qts(
        const Rcpp::DataFrame &qts,
        const bool disable_normalization = false
);

//' QTS Normalization
//'
//' This function ensures that all quaternions in the time series are unit
//' quaternions.
//'
//' @param qts A quaternion time series stored as a \code{\link[tibble]{tibble}}
//'   with columns `time`, `w`, `x`, `y` and `z`.
//'
//' @return A quaternion time series stored as a \code{\link[tibble]{tibble}}
//'   with columns `time`, `w`, `x`, `y` and `z` in which quaternions are unit
//'   quaternions.
//'
//' @export
//' @examples
//' # TO DO
// [[Rcpp::export]]
Rcpp::DataFrame normalize_qts(const Rcpp::DataFrame &qts);

// [[Rcpp::export]]
Rcpp::DataFrame derivative_qts_impl(const Rcpp::DataFrame &qts);

// [[Rcpp::export]]
Rcpp::DataFrame log_qts_impl(const Rcpp::DataFrame &qts);

// [[Rcpp::export]]
Rcpp::DataFrame exp_qts_impl(const Rcpp::DataFrame &qts);

//' QTS Centering and Standardization
//'
//' This function operates a centring of the QTS around the geometric mean of
//' its quaternions. This is effectively achieved by premultiplying each
//' quaternion by the inverse of their geometric mean.
//'
//' @param qts A quaternion time series stored as a \code{\link[tibble]{tibble}}
//'   with columns `time`, `w`, `x`, `y` and `z`.
//' @param standardize A boolean specifying whether to standardize the QTS in
//'   addition to centering it. Defaults to `FALSE`.
//'
//' @return A quaternion time series stored as a \code{\link[tibble]{tibble}}
//'   with columns `time`, `w`, `x`, `y` and `z` in which quaternions have been
//'   centered around their geometric mean.
//'
//' @export
//' @examples
//' # TO DO
// [[Rcpp::export]]
Rcpp::List centring_qts(const Rcpp::DataFrame &qts, const bool standardize = false);

#endif /* SQUATQTSCLASS_H */
