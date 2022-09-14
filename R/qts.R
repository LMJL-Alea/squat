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
#'   attitude at current time point.
#'
#' @export
#' @examples
#' derivative_qts(vespa$igp[[1]])
derivative_qts <- function(qts) {
  qts <- derivative_qts_impl(qts)
  qts[-1, ]
}

#' QTS Random Sampling
#'
#' This function adds uncorrelated Gaussian noise to the logarithm QTS using an
#' exponential covariance function.
#'
#' See \code{\link[roahd]{exp_cov_function}} for details about the roles of
#' `alpha` and `beta` in the definition of the covariance operator.
#'
#' @param n An integer specifying how many QTS should be generated.
#' @param mean A \code{\link[tibble]{tibble}} specifying the mean QTS.
#' @param alpha A positive scalar specifying the variance of each component of
#'   the log-QTS. Defaults to `0.01`.
#' @param beta A positive scalar specifying the exponential weight. Defaults to
#'   `0.001`.
#'
#' @return A list of `n` QTS with added noise as specified by parameters `alpha`
#'   and `beta`.
#' @export
#'
#' @examples
#' rnorm_qts(1, vespa$igp[[1]])
rnorm_qts <- function(n, mean, alpha = 0.01, beta = 0.001) {
  mean <- log_qts(mean)
  time_grid <- mean$time
  C1 <- roahd::exp_cov_function(time_grid, alpha = alpha, beta = beta)
  C2 <- roahd::exp_cov_function(time_grid, alpha = alpha, beta = beta)
  C3 <- roahd::exp_cov_function(time_grid, alpha = alpha, beta = beta)
  centerline <- rbind(mean$x, mean$y, mean$z)
  CholC1 <- chol(C1)
  CholC2 <- chol(C2)
  CholC3 <- chol(C3)
  mean <- roahd::generate_gauss_mfdata(
    N = n,
    L = 3,
    centerline = centerline,
    correlations = rep(0, 3),
    listCholCov = list(CholC1, CholC2, CholC3)
  )
  qts_list <- purrr::map(1:n, ~ {
    tibble(
      time = time_grid,
      w = 0,
      x = mean[[1]][.x, ],
      y = mean[[2]][.x, ],
      z = mean[[3]][.x, ]
    )
  })
  purrr::map(qts_list, exp_qts)
}

#' QTS Sample Centering and Standardization
#'
#' @param qts_list A list of quaternion time series stored as
#'   \code{\link[tibble]{tibble}}s with columns `time`, `w`, `x`, `y` and `z`.
#' @param center A boolean specifying whether to center the sample of QTS or
#'   not. If set to `FALSE`, the original sample is returned, meaning that no
#'   standardization is performed regardless of whether argument `standardize`
#'   was set to `TRUE` or not. Defaults to `TRUE`.
#' @param standardize A boolean specifying whether to standardize the sample of
#'   QTS once they have been centered. Defaults to `TRUE`.
#' @param by_row A boolean specifying whether the QTS scaling should happen for
#'   each data point (`by_row = TRUE`) or for each time point (`by_row =
#'   FALSE`). Defaults to `FALSE`.
#' @param keep_summary_stats A boolean specifying whether the mean and standard
#'   deviation used for standardizing the data. Defaults to `FALSE` in which
#'   case only the list of properly rescaled QTS is returned.
#'
#' @return A list of properly rescaled QTS when `keep_summary_stats = FALSE`.
#'   Otherwise a list with three components:
#' - `qts_list`: a list of properly rescaled QTS;
#' - `mean`: a numeric vector with the quaternion Fréchet mean;
#' - `sd`: a numeric vector with the quaternion Fréchet standard deviation.
#'
#' @export
#'
#' @examples
#' qts_list <- scale_qts(vespa$igp)
#' qts_list[[1]]
scale_qts <- function(qts_list,
                      center = TRUE,
                      standardize = TRUE,
                      by_row = FALSE,
                      keep_summary_stats = FALSE) {
  if (!center) {
    if (!keep_summary_stats) return(qts_list)
    return(list(
      qts_list = qts_list,
      mean_values = NA,
      sd_values = NA
    ))
  }

  if (!by_row) {
    qts_list <- qts_list |>
      purrr::map(purrr::array_tree, margin = 1) |>
      purrr::transpose() |>
      purrr::map(purrr::reduce, rbind) |>
      purrr::map(tibble::as_tibble)
  }

  std_data <- purrr::map(qts_list, centring_qts, standardize = standardize)
  qts_list <- purrr::map(std_data, "qts")

  if (!by_row) {
    qts_list <- qts_list |>
      purrr::map(purrr::array_tree, margin = 1) |>
      purrr::transpose() |>
      purrr::map(purrr::reduce, rbind) |>
      purrr::map(tibble::as_tibble)
  }

  if (!keep_summary_stats) return(qts_list)

  list(
    qts_list = qts_list,
    mean_values = purrr::map(std_data, "mean"),
    sd_values = purrr::map_dbl(std_data, "sd")
  )
}

#' QTS Straightening
#'
#' This function straightens a QTS so that the last point equals the first
#' point.
#'
#' @param qts A \code{\link[tibble]{tibble}} specifying a QTS.
#'
#' @return A \code{\link[tibble]{tibble}} storing the straightened QTS.
#' @export
#'
#' @examples
#' straighten_qts(vespa$igp[[1]])
straighten_qts <- function(qts) {
  P <- nrow(qts)
  time_values <- qts$time
  t1 <- time_values[1]
  tP <- time_values[P]
  qts <- log_qts(qts)
  for (i in 3:5) {
    y1 <- qts[[i]][1]
    yP <- qts[[i]][P]
    a <- (yP - y1) / (tP - t1)
    qts[[i]] <- qts[[i]] - a * (time_values - t1)
  }
  qts$w <- rep(0, P)
  exp_qts(qts)
}
