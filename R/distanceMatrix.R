#' Distance Matrix for Quaternion Time Series Samples
#'
#' @param qts_list A list of \code{\link[tibble]{tibble}}s storing a sample of
#'   quaternion time series.
#' @param normalize_distance A boolean specifying whether to compute normalized
#'   distance between QTS. Please note that not all step patterns are
#'   normalizable. Defaults to `FALSE`.
#' @param labels A character vector specifying labels for each QTS. Defaults to
#'   `NULL` which uses row numbers as labels.
#' @param quotient_space A boolean specifying whether to make the distance
#'   invariant to the coordinate system into which rotations are expressed.
#'   Defaults to `FALSE` for backward compatibility and faster computation.
#' @inheritParams DTW
#'
#' @return A \code{\link[stats]{dist}} object storing the distance matrix
#'   between QTS in a sample via DTW.
#' @export
#'
#' @examples
#' # TO DO
distDTW <- function(qts_list,
                    normalize_distance = TRUE,
                    labels = NULL,
                    resample = TRUE,
                    disable_normalization = FALSE,
                    step_pattern = dtw::symmetric2,
                    quotient_space = FALSE) {
  if (normalize_distance && is.na(attr(step_pattern, "norm")))
    stop("The provided step pattern is not normalizable.")

  if (!disable_normalization) {
    qts_list <- purrr::map(qts_list, normalize_qts)
  }

  if (resample) {
    qts_list <- purrr::map(qts_list, resample_qts, disable_normalization = TRUE)
  }

  n <- length(qts_list)
  if (is.null(labels))
    labels <- 1:n

  indices <- linear_index(n)

  if (quotient_space) {
    d <- furrr::future_map_dbl(indices$k, ~ {
      i <- indices$i[.x]
      j <- indices$j[.x]
      dtw_data <- DTWi(
        qts1 = qts_list[[i]],
        qts2 = qts_list[[j]],
        resample = FALSE,
        disable_normalization = TRUE,
        distance_only = TRUE,
        step_pattern = step_pattern
      )
      if (normalize_distance)
        dtw_data$normalizedDistance
      else
        dtw_data$distance
    }, .options = furrr::furrr_options(seed = TRUE))
  } else {
    d <- furrr::future_map_dbl(indices$k, ~ {
      i <- indices$i[.x]
      j <- indices$j[.x]
      dtw_data <- DTW(
        qts1 = qts_list[[i]],
        qts2 = qts_list[[j]],
        resample = FALSE,
        disable_normalization = TRUE,
        distance_only = TRUE,
        step_pattern = step_pattern
      )
      if (normalize_distance)
        dtw_data$normalizedDistance
      else
        dtw_data$distance
    }, .options = furrr::furrr_options(seed = TRUE))
  }

  attributes(d) <- NULL
  attr(d, "Labels") <- labels
  attr(d, "Size") <- n
  attr(d, "call") <- match.call()
  class(d) <- "dist"
  d
}

linear_index <- function(n) {
  res <- purrr::cross_df(
    .l = list(j = 1:n, i = 1:n),
    .filter = ~ .x <= .y
  )
  res$k <- n * (res$i - 1) - res$i * (res$i - 1) / 2 + res$j - res$i
  res
}
