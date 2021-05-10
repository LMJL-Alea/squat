#' QTS Derivative
#'
#' This function computes the first derivative of a quaternion time series with
#' respect to time.
#'
#' @param qts A quaternion time series stored as a \code{\link[tibble]{tibble}}
#'   with columns `time`, `w`, `x`, `y` and `z`.
#'
#' @return A quaternion time series stored as a \code{\link[tibble]{tibble}}
#'   with columns `time`, `w`, `x`, `y` and `z` in which quaternions measure
#'   the rotation to be applied to transform attitude at previous time point to
#'   attitude at the current time point.
#'
#' @export
#' @examples
#' TO DO
derivative_qts <- function(qts) {
  qts <- derivative_qts_impl(qts)
  qts[-1, ]
}
