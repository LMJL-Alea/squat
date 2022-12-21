#' Distance Matrix for Quaternion Time Series Samples
#'
#' @param qts_list An object of class [qts_sample].
#' @param normalize_distance A boolean specifying whether to compute normalized
#'   distance between QTS. Please note that not all step patterns are
#'   normalizable. Defaults to `FALSE`.
#' @param labels A character vector specifying labels for each QTS. Defaults to
#'   `NULL` which uses row numbers as labels.
#' @inheritParams DTW
#'
#' @return A [stats::dist] object storing the distance matrix between QTS in a
#'   sample via DTW.
#' @export
#'
#' @examples
#' D <- distDTW(vespa64$igp)
distDTW <- function(qts_list,
                    normalize_distance = TRUE,
                    labels = NULL,
                    resample = TRUE,
                    disable_normalization = FALSE,
                    step_pattern = dtw::symmetric2) {
  if (!is_qts_sample(qts_list))
    cli::cli_abort("The input argument {.arg qts_list} should be of class {.cls qts_sample}. You can try {.fn as_qts_sample()}.")

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

  .pairwise_distances <- function(linear_indices) {
    pb <- progressr::progressor(along = linear_indices)
    furrr::future_map_dbl(linear_indices, ~ {
      pb()
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

  d <- .pairwise_distances(indices$k)

  attributes(d) <- NULL
  attr(d, "Labels") <- labels
  attr(d, "Size") <- n
  attr(d, "call") <- match.call()
  class(d) <- "dist"
  d
}

linear_index <- function(n) {
  res <- tidyr::expand_grid(i = 1:n, j = 1:n)
  res <- subset(res, res$j > res$i)
  res$k <- n * (res$i - 1) - res$i * (res$i - 1) / 2 + res$j - res$i
  res
}
