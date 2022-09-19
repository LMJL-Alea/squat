#ifndef SQUATQTSCLASS_H
#define SQUATQTSCLASS_H

#include <Rcpp.h>

//' QTS Transformation To Distance Time Series
//'
//' This function computes a real-valued time series reporting the pointwise
//' geodesic distance between the two input QTS at each time point.
//'
//' The function currently expects that the two input QTS are evaluated on the
//' same time grid and does not check this assumption.
//'
//' @param first_qts A quaternion time series stored as a
//'   \code{\link[tibble]{tibble}} with columns `time`, `w`, `x`, `y` and `z`.
//' @param second_qts A quaternion time series stored as a
//'   \code{\link[tibble]{tibble}} with columns `time`, `w`, `x`, `y` and `z`.
//'
//' @return A time series stored as a \code{\link[tibble]{tibble}} with columns
//'   `time` and `distance` in which `distance` measures the angular distance
//'   between the quaternions of both input QTS at a given time point.
//'
//' @export
//' @examples
//' # TO DO
// [[Rcpp::export]]
Rcpp::DataFrame qts2distance(
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

#endif /* SQUATQTSCLASS_H */
