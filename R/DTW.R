#' Dynamic Time Warping for Quaternion Time Series
#'
#' This function evaluates the Dynamic Time Warping (DTW) distance between two
#' quaternion time series (QTS).
#'
#' If no evaluation grid is provided, the function assumes that the two input
#' QTS are evaluated on the same grid.
#'
#' @param qts1 A \code{\link[tibble]{tibble}} storing the 1st quaternion time
#'   series.
#' @param qts2 A \code{\link[tibble]{tibble}} storing the 2nd quaternion time
#'   series.
#' @param resample A boolean specifying whether the QTS should be uniformly
#'   resampled on their domain before computing distances. Defaults to `TRUE`.
#' @param disable_normalization A boolean specifying whether quaternion
#'   normalization should be disabled. Defaults to `FALSE` which ensures that we
#'   always deal with unit quaternions.
#' @param distance_only A boolean specifying whether to only compute distance
#'   (no backtrack, faster). Defaults to `FALSE`.
#' @param step_pattern A \code{\link[dtw]{stepPattern}} specifying the local
#'   constraints on the warping path. Defaults to \code{\link[dtw]{symmetric2}}
#'   which uses symmetric and normalizable warping paths with no local slope
#'   constraints. See \code{\link[dtw]{stepPattern}} for more information.
#'
#' @return An object of class \code{\link[dtw]{dtw}} storing the dynamic time
#'   warping results.
#' @export
#'
#' @examples
#' TO DO
DTW <- function(qts1, qts2,
                resample = TRUE,
                disable_normalization = FALSE,
                distance_only = FALSE,
                step_pattern = dtw::symmetric2) {
  if (!disable_normalization) {
    qts1 <- normalize_qts(qts1)
    qts2 <- normalize_qts(qts2)
  }
  if (resample) {
    qts1 <- resample_qts(qts1, disable_normalization = TRUE)
    qts2 <- resample_qts(qts2, disable_normalization = TRUE)
  }
  M <- GetCostMatrix(qts1, qts2, disable_normalization = TRUE)
  dtw::dtw(M, distance.only = distance_only, step.pattern = step_pattern)
}
