#' QTS Sample Class
#'
#' A collection of functions that implements the QTS sample class. It currently
#' provides the `as_qts_sample()` function for QTS sample coercion of lists of
#' \code{\link{qts}} objects and the `is_qts_sample()` function for checking if
#' an object is a QTS sample.
#'
#' A QTS sample is a collection of quaternion time series (QTS), each of which
#' is stored as a \code{\link[tibble]{tibble}} with 5 columns:
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
#' @param x A list of \code{\link[tibble]{tibble}}s, each of which with columns
#'   `time`, `w`, `x`, `y` and `z`.
#'
#' @return An object of class \code{\link{qts_sample}}.
#' @name qts_sample
#'
#' @examples
#' x <- vespa$igp
#' y <- as_qts_sample(x)
#' is_qts_sample(x)
#' is_qts_sample(y)
NULL

#' @export
#' @rdname qts_sample
as_qts_sample <- function(x) {
  if (!is.list(x))
    cli::cli_abort("The input {.arg x} should be a list.")
  for (qts in x) {
    if (!is_qts(qts))
      cli::cli_abort("All elements of the input list should be of class {.cls qts}.")
  }
  class(x) <- c("qts_sample", class(x))
  x
}

#' @export
#' @rdname qts_sample
is_qts_sample <- function(x) {
  "qts_sample" %in% class(x)
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
#' rnorm_qts(1, vespa64$igp[[1]])
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
  out <- purrr::map(qts_list, exp_qts)
  as_qts_sample(out)
}

#' QTS Sample Centering and Standardization
#'
#' @param x An object of class \code{\link{qts_sample}} representing a sample of
#'   observed QTS.
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
#'   deviation used for standardizing the data should be stored in the output
#'   object. Defaults to `FALSE` in which case only the list of properly
#'   rescaled QTS is returned.
#'
#' @return A list of properly rescaled QTS stored as an object of class
#'   \code{\link{qts_sample}} when `keep_summary_stats = FALSE`.
#'   Otherwise a list with three components:
#' - `rescaled_sample`: a list of properly rescaled QTS stored as an object of
#' class \code{\link{qts_sample}};
#' - `mean`: a numeric vector with the quaternion Fréchet mean;
#' - `sd`: a numeric vector with the quaternion Fréchet standard deviation.
#'
#' @export
#'
#' @examples
#' x <- scale_qts(vespa64$igp)
#' x[[1]]
scale_qts <- function(x,
                      center = TRUE,
                      standardize = TRUE,
                      by_row = FALSE,
                      keep_summary_stats = FALSE) {
  if (!is_qts_sample(x))
    cli::cli_abort("The input argument {.arg x} should be of class {.cls qts_sample}.")

  if (!center) {
    if (!keep_summary_stats) return(x)
    return(list(
      rescaled_sample = x,
      mean_values = NA,
      sd_values = NA
    ))
  }

  if (!by_row) {
    x <- x |>
      purrr::map(purrr::array_tree, margin = 1) |>
      purrr::transpose() |>
      purrr::map(purrr::reduce, rbind) |>
      purrr::map(tibble::as_tibble) |>
      purrr::map(as_qts)
  }

  std_data <- purrr::map(x, centring_qts,
                         standardize = standardize,
                         keep_summary_stats = TRUE)
  x <- purrr::map(std_data, "qts")

  if (!by_row) {
    x <- x |>
      purrr::map(purrr::array_tree, margin = 1) |>
      purrr::transpose() |>
      purrr::map(purrr::reduce, rbind) |>
      purrr::map(tibble::as_tibble) |>
      purrr::map(as_qts)
  }

  if (!keep_summary_stats) return(as_qts_sample(x))

  list(
    rescaled_sample = as_qts_sample(x),
    mean_values = purrr::map(std_data, "mean"),
    sd_values = purrr::map_dbl(std_data, "sd")
  )
}

#' QTS Geometric Mean
#'
#' This function computes the pointwise geometric mean of a QTS sample.
#'
#' @param x An object of class \code{\link{qts_sample}}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An object of class \code{\link{qts}} in which quaternions are the
#'   pointwise geometric mean of the input QTS sample.
#'
#' @export
#' @examples
#' mean(vespa64$igp)
mean.qts_sample <- function(x, ...) {
  if (!is_qts_sample(x)) x <- as_qts_sample(x)
  mean_qts_impl(x)
}

#' QTS Geometric Median
#'
#' This function computes the pointwise geometric median of a QTS sample.
#'
#' @param x An object of class \code{\link{qts_sample}}.
#' @param na.rm A logical value indicating whether NA values should be stripped
#'   before the computation proceeds.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An object of class \code{\link{qts}} in which quaternions are the
#'   pointwise geometric median of the input QTS sample.
#'
#' @importFrom stats median
#' @export
#' @examples
#' median(vespa64$igp)
median.qts_sample <- function(x, na.rm = FALSE, ...) {
  if (!is_qts_sample(x)) x <- as_qts_sample(x)
  median_qts_impl(x)
}
