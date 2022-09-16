#' QTS Class
#'
#' A collection of functions that implements the QTS class. It currently
#' provides the `as_qts()` function for QTS coercion of tibbles and the
#' `is_qts()` function for checking is an object is a QTS.
#'
#' A quaternion time series (QTS) is stored as a \code{\link[tibble]{tibble}}
#' with 5 columns:
#' - `time`: A first column specifying the time points at which quaternions were
#' collected;
#' - `w`: A second column specifying the first coordinate of the collected
#' quaternions;
#' - `x`: A third column specifying the second coordinate of the collected
#' quaternions;
#' - `y`: A fourth column specifying the third coordinate of the collected
#' quaternions;
#' - `z`: A fifth column specifying the fourth coordinate of the collected
#' quaternions.
#'
#' @param x A \code{\link[tibble]{tibble}} with columns `time`, `w`, `x`, `y`
#'   and `z`.
#'
#' @return An object of class \code{\link{qts}}.
#' @name qts
#'
#' @examples
#' qts1 <- vespa$igp[[1]]
#' qts2 <- as_qts(qts1)
#' is_qts(qts1)
#' is_qts(qts2)
NULL

#' @export
#' @rdname qts
as_qts <- function(x) {
  if (is_qts(x)) return(x)
  if (!tibble::is_tibble(x))
    cli::cli_abort("The input object should be of class {.cls tbl}.")
  if (!all(names(x) == c("time", "w", "x", "y", "z")))
    cli::cli_abort("The input tibble should have exactly the 5 following columns in that order: {.code time}, {.code w}, {.code x}, {.code y} and {.code z}.")
  class(x) <- c("qts", class(x))
  x
}

#' @export
#' @rdname qts
is_qts <- function(x) {
  "qts" %in% class(x)
}

#' QTS Derivative
#'
#' This function computes the first derivative of a quaternion time series with
#' respect to time.
#'
#' @param qts An object of class \code{\link{qts}}.
#'
#' @return An object of class \code{\link{qts}} in which quaternions measure
#'   the rotation to be applied to transform attitude at previous time point to
#'   attitude at current time point.
#'
#' @export
#' @examples
#' derivative_qts(vespa$igp[[1]])
derivative_qts <- function(qts) {
  if (!is_qts(qts))
    cli::cli_abort("The input object should be of class {.cls qts}.")
  qts <- derivative_qts_impl(qts)
  as_qts(qts[-1, ])
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
#' @param mean_qts An object of class \code{\link{qts}} specifying the mean QTS.
#' @param alpha A positive scalar specifying the variance of each component of
#'   the log-QTS. Defaults to `0.01`.
#' @param beta A positive scalar specifying the exponential weight. Defaults to
#'   `0.001`.
#'
#' @return A list of `n` objects of class \code{\link{qts}} with added noise as
#'   specified by parameters `alpha` and `beta`.
#' @export
#'
#' @examples
#' rnorm_qts(1, vespa$igp[[1]])
rnorm_qts <- function(n, mean_qts, alpha = 0.01, beta = 0.001) {
  if (!is_qts(mean_qts))
    cli::cli_abort("The {.arg mean_qts} parameter should be of class {.cls qts}.")
  mean_qts <- log_qts(mean_qts)
  time_grid <- mean_qts$time
  C1 <- roahd::exp_cov_function(time_grid, alpha = alpha, beta = beta)
  C2 <- roahd::exp_cov_function(time_grid, alpha = alpha, beta = beta)
  C3 <- roahd::exp_cov_function(time_grid, alpha = alpha, beta = beta)
  centerline <- rbind(mean_qts$x, mean_qts$y, mean_qts$z)
  CholC1 <- chol(C1)
  CholC2 <- chol(C2)
  CholC3 <- chol(C3)
  mean_qts <- roahd::generate_gauss_mfdata(
    N = n,
    L = 3,
    centerline = centerline,
    correlations = rep(0, 3),
    listCholCov = list(CholC1, CholC2, CholC3)
  )
  qts_list <- purrr::map(1:n, ~ {
    as_qts(tibble(
      time = time_grid,
      w = 0,
      x = mean_qts[[1]][.x, ],
      y = mean_qts[[2]][.x, ],
      z = mean_qts[[3]][.x, ]
    ))
  })
  purrr::map(qts_list, exp_qts)
}

#' QTS Sample Centering and Standardization
#'
#' @param qts_list A list of objects of class \code{\link{qts}} representing a
#'   sample of QTS observations.
#' @param center A boolean specifying whether to center the sample of QTS. If
#'   set to `FALSE`, the original sample is returned, meaning that no
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
#' @return A list of properly rescaled QTS stored as objects of class
#'   \code{\link{qts}} when `keep_summary_stats = FALSE`.
#'   Otherwise a list with three components:
#' - `qts_list`: a list of properly rescaled QTS stored as objects of class
#' \code{\link{qts}};
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
  for (qts in qts_list) {
    if (!is_qts(qts))
      cli::cli_abort("All objects in the input sample {.arg qts_list} should be of class {.cls qts}.")
  }

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
    qts_list = purrr::map(qts_list, as_qts),
    mean_values = purrr::map(std_data, "mean"),
    sd_values = purrr::map_dbl(std_data, "sd")
  )
}

#' QTS Straightening
#'
#' This function straightens a QTS so that the last point equals the first
#' point.
#'
#' @param qts An object of class \code{\link{qts}}.
#'
#' @return An object of class \code{\link{qts}} storing the straightened QTS.
#' @export
#'
#' @examples
#' straighten_qts(vespa$igp[[1]])
straighten_qts <- function(qts) {
  if (!is_qts(qts))
    cli::cli_abort("The input object should be of class {.cls qts}.")
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
