#ifndef SQUATQTSTRANSFORMATIONS_H
#define SQUATQTSTRANSFORMATIONS_H

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::DataFrame qts2dts_impl(
    const Rcpp::DataFrame &first_qts,
    const Rcpp::DataFrame &second_qts
);

//' QTS Transformation To Norm Time Series
//'
//' This function computes a univariate time series representing the norm of the
//' quaternions.
//'
//' @param qts A quaternion time series stored as a \code{\link[tibble]{tibble}}
//'   with columns `time`, `w`, `x`, `y` and `z`.
//' @param disable_normalization A boolean specifying whether quaternion
//'   normalization should be disabled. Defaults to `FALSE`.
//'
//' @return A time series stored as a \code{\link[tibble]{tibble}} with columns
//'   `time` and `norm` in which `norm` measures the angular distance between
//'   the current quaternion and the identity.
//'
//' @export
//' @examples
//' # TO DO
// [[Rcpp::export]]
Rcpp::DataFrame qts2norm(
    const Rcpp::DataFrame &qts,
    const bool disable_normalization = false
);

//' QTS Transformation To Angle Time Series
//'
//' This function computes a univariate time series representing the angle
//' between the first and other attitudes.
//'
//' @param qts A quaternion time series stored as a \code{\link[tibble]{tibble}}
//'   with columns `time`, `w`, `x`, `y` and `z`.
//' @param disable_normalization A boolean specifying whether quaternion
//'   normalization should be disabled. Defaults to `FALSE`.
//'
//' @return A time series stored as a \code{\link[tibble]{tibble}} with columns
//'   `time` and `angle` in which `angle` measures the angle between the current
//'   rotation and the first one.
//'
//' @export
//' @examples
//' # TO DO
// [[Rcpp::export]]
Rcpp::DataFrame qts2angle(
    const Rcpp::DataFrame &qts,
    const bool disable_normalization = false
);

//' QTS Transformation to Angular Velocity Time Series
//'
//' This function projects a quaternion time series into the space of angular
//' velocities.
//'
//' @param qts A QTS stored as a \code{\link[tibble]{tibble}}s with columns
//'   `time`, `w`, `x`, `y` and `z`.
//' @param fixed_frame A string specifying the fixed frame with respect to which
//'   coordinates of the angular velocity should be computed. Choices are
//'   `"global"` or `"body"`. Defaults to `"global"`.
//'
//' @return A time series stored as a \code{\link[tibble]{tibble}} with columns
//'   `time`, `x`, `y` and `z` containing the angular velocity at each time
//'   point.
//'
//' @export
//' @examples
//' # TO DO
// [[Rcpp::export]]
Rcpp::DataFrame qts2avts(
    const Rcpp::DataFrame &qts,
    const Rcpp::String &fixed_frame = "global"
);

//' QTS Transformation from Angular Velocity Time Series
//'
//' This function projects back an angular velocity time series into the space
//' of quaternions.
//'
//' @param avts An angular velocity time series stored as a
//'   \code{\link[tibble]{tibble}}s with columns `time`, `x`, `y` and `z`.
//' @param init_t A positive scalar specifying the initial time point.
//' @param init_q A length-4 numeric vector speicyfing the quaternion at the
//'   initial time point.
//'
//' @return A quaternion time series stored as a \code{\link[tibble]{tibble}}
//'   with columns `time`, `w`, `x`, `y` and `z`.
//'
//' @export
//' @examples
//' # TO DO
// [[Rcpp::export]]
Rcpp::DataFrame avts2qts(
    const Rcpp::DataFrame &avts,
    const double init_t,
    const Rcpp::NumericVector init_q
);

#endif /* SQUATQTSTRANSFORMATIONS_H */
