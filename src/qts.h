#ifndef QTS_H
#define QTS_H

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
//' TO DO
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
//' TO DO
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
//' TO DO
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
//' TO DO
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
//' TO DO
// [[Rcpp::export]]
Rcpp::DataFrame normalize_qts(const Rcpp::DataFrame &qts);

// [[Rcpp::export]]
Rcpp::DataFrame derivative_qts_impl(const Rcpp::DataFrame &qts);

//' QTS Logarithm
//'
//' This function computes the logarithm of a quaternion time series as the time
//' series of the quaternion logarithms.
//'
//' @param qts A quaternion time series stored as a \code{\link[tibble]{tibble}}
//'   with columns `time`, `w`, `x`, `y` and `z`.
//'
//' @return A quaternion time series stored as a \code{\link[tibble]{tibble}}
//'   with columns `time`, `w`, `x`, `y` and `z` in which quaternions have been
//'   replaced by their logarithm.
//'
//' @export
//' @examples
//' TO DO
// [[Rcpp::export]]
Rcpp::DataFrame log_qts(const Rcpp::DataFrame &qts);

//' QTS Exponential
//'
//' This function computes the exponential of a quaternion time series as the
//' time series of the quaternion exponentials.
//'
//' @param qts A quaternion time series stored as a \code{\link[tibble]{tibble}}
//'   with columns `time`, `w`, `x`, `y` and `z`.
//'
//' @return A quaternion time series stored as a \code{\link[tibble]{tibble}}
//'   with columns `time`, `w`, `x`, `y` and `z` in which quaternions have been
//'   replaced by their exponential.
//'
//' @export
//' @examples
//' TO DO
// [[Rcpp::export]]
Rcpp::DataFrame exp_qts(const Rcpp::DataFrame &qts);

//' @export
// [[Rcpp::export]]
Rcpp::DataFrame centring_qts(const Rcpp::DataFrame &qts);

#endif /* QTS_H */
