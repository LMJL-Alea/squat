#' Distance Matrix for Quaternion Time Series Samples
#'
#' @param qts_list A list of \code{\link[tibble]{tibble}}s storing a sample of quaternion time series.
#' @param normalize_distance A boolean specifying whether to compute normalized
#'   distance between QTS. Please note that not all step patterns are
#'   normalizable. Defaults to `FALSE`.
#' @param labels A character vector specifying labels for each QTS. Defaults to `NULL` which uses row numbers as labels.
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
                    step_pattern = dtw::symmetric2) {
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

  d <- numeric(n * (n - 1)/2)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      dtw_data <- DTW(
        qts1 = qts_list[[i]],
        qts2 = qts_list[[j]],
        resample = FALSE,
        disable_normalization = TRUE,
        distance_only = TRUE,
        step_pattern = step_pattern
      )
      d[n * (i - 1) - i * (i - 1)/2 + j - i] <- if (normalize_distance)
        dtw_data$normalizedDistance
      else
        dtw_data$distance
    }
  }
  attributes(d) <- NULL
  attr(d, "Labels") <- labels
  attr(d, "Size") <- n
  attr(d, "call") <- match.call()
  class(d) <- "dist"
  d
}
