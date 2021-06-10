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
#' # TO DO
derivative_qts <- function(qts) {
  qts <- derivative_qts_impl(qts)
  qts[-1, ]
}

#' QTS Noise Generator
#'
#' This function adds uncorrelated Gaussian noise to the logarithm QTS using an
#' exponential covariance function.
#'
#' See \code{\link[roahd]{exp_cov_function}} for details about the roles of
#' `alpha` and `beta` in the definition of the covariance operator.
#'
#' @param qts A \code{\link[tibble]{tibble}} specifying the input QTS.
#' @param n An integer specifying how many noisy QTS should be generated.
#'   Defaults to `1L`.
#' @param alpha A positive scalar specifying the variance of each component of
#'   the log-QTS. Defaults to `0.2`.
#' @param beta A positive scalar specifying the exponential weight. Defaults to
#'   `0.5`.
#'
#' @return A list of `n` QTS with added noise as specified by parameters `alpha`
#'   and `beta`.
#' @export
#'
#' @examples
#' qts <- tibble::tibble(
#'   time = seq(0, 1, len = 101),
#'   w = 0,
#'   x = cos(2 * pi * time),
#'   y = time^2,
#'   z = sin(time)
#' )
#' qts <- exp_qts(qts)
#' qts_n <- add_noise(qts, n = 10)
add_noise <- function(qts, n = 1, alpha = 0.1, beta = 0.5) {
  qts <- log_qts(qts)
  time_grid <- qts$time
  C1 <- roahd::exp_cov_function(time_grid, alpha = alpha, beta = beta)
  C2 <- roahd::exp_cov_function(time_grid, alpha = alpha, beta = beta)
  C3 <- roahd::exp_cov_function(time_grid, alpha = alpha, beta = beta)
  centerline <- rbind(qts$x, qts$y, qts$z)
  CholC1 <- chol(C1)
  CholC2 <- chol(C2)
  CholC3 <- chol(C3)
  qts <- roahd::generate_gauss_mfdata(
    N = n,
    L = 3,
    centerline = centerline,
    correlations = rep(0, 3),
    listCholCov = list(CholC1, CholC2, CholC3)
  )
  qts_list <- lapply(1:n, function(.x) {
    tibble(
      time = time_grid,
      w = 0,
      x = qts[[1]][.x, ],
      y = qts[[2]][.x, ],
      z = qts[[3]][.x, ]
    )
  })
  lapply(qts_list, exp_qts)
}
