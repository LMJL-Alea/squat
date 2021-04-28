#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <Rcpp.h>

// [[Rcpp::interfaces(r, cpp)]]

//' QTS Resampling
//'
//' This function performs uniform resampling using SLERP.
//'
//' @param qts A quaternion time series stored as a \code{\link[tibble]{tibble}}
//'   with columns `time`, `w`, `x`, `y` and `z`.
//' @param tmin A numeric value specifying the lower bound of the time interval
//'   over which uniform resampling should take place. It must satisfy `tmin >=
//'   min(qts$time)`. Defaults to `NA` in which case it is set to
//'   `min(qts$time)`.
//' @param tmax A numeric value specifying the upper bound of the time interval
//'   over which uniform resampling should take place. It must satisfy `tmax <=
//'   max(qts$time)`. Defaults to `NA` in which case it is set to
//'   `max(qts$time)`.
//' @param nout An integer specifying the size of the uniform grid for time
//'   resampling. Defaults to `0L` in which case it uses the same grid size as
//'   the input QTS.
//' @param disable_normalization A boolean specifying whether quaternion
//'   normalization should be disabled. Defaults to `FALSE` in which case the
//'   function makes sure that quaternions are normalized prior to performing
//'   SLERP interpolation.
//'
//' @return A quaternion time series stored as a \code{\link[tibble]{tibble}}
//'   with columns `time`, `w`, `x`, `y` and `z` in which quaternions are
//'   uniformly sampled in the range `[tmin, tmax]`.
//'
//' @export
//' @examples
//' TO DO
// [[Rcpp::export]]
Rcpp::DataFrame resample_qts(
    const Rcpp::DataFrame &qts,
    double tmin = NA_REAL,
    double tmax = NA_REAL,
    const unsigned int nout = 0,
    const bool disable_normalization = false
);

//' @export
// [[Rcpp::export]]
Rcpp::DataFrame smooth_qts(
    const Rcpp::DataFrame &qts,
    const double alpha = 0.5
);

#endif /* INTERPOLATION_H */
