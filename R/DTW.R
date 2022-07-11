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
#' # TO DO
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

DTWi <- function(qts1, qts2,
                 resample = TRUE,
                 disable_normalization = FALSE,
                 distance_only = FALSE,
                 step_pattern = dtw::symmetric2) {
  opt <- nloptr::directL(
    fn = cost,
    lower = rep(0, 3),
    upper = c(2 * pi, pi, 2 * pi),
    qts1 = qts1, qts2 = qts2,
    resample = resample,
    disable_normalization = disable_normalization,
    distance_only = TRUE,
    step_pattern = step_pattern
  )
  if (anyNA(opt$par)) {
    cli::cli_alert_danger("Optimization failed")
    return(DTW(qts1, qts2,
        resample = TRUE,
        disable_normalization = FALSE,
        distance_only = FALSE,
        step_pattern = dtw::symmetric2))
  }
  for (i in 1:3) {
    if (opt$par[1] < 0)
      opt$par[1] <- .Machine$double.eps^0.5
  }
  if (opt$par[1] > 2*pi)
    opt$par[1] <- 2*pi - .Machine$double.eps^0.5
  if (opt$par[2] > pi)
    opt$par[2] <- pi - .Machine$double.eps^0.5
  if (opt$par[3] > 2*pi)
    opt$par[3] <- 2*pi - .Machine$double.eps^0.5
  opt <- nloptr::bobyqa(
    x0 = opt$par,
    fn = cost,
    lower = rep(0, 3),
    upper = c(2 * pi, pi, 2 * pi),
    qts1 = qts1, qts2 = qts2,
    resample = resample,
    disable_normalization = disable_normalization,
    distance_only = TRUE,
    step_pattern = step_pattern
  )
  cost_impl(opt$par,
            qts1, qts2,
            resample,
            disable_normalization,
            distance_only,
            step_pattern)
}

cost <- function(q0,
                 qts1, qts2,
                 resample,
                 disable_normalization,
                 distance_only,
                 step_pattern) {
  cost_impl(
    q0,
    qts1, qts2,
    resample,
    disable_normalization,
    distance_only,
    step_pattern
  )$distance
}

cost_impl <- function(q0,
                      qts1, qts2,
                      resample,
                      disable_normalization,
                      distance_only,
                      step_pattern) {
  omega <- q0[1]
  theta <- q0[2]
  phi <- q0[3]
  q0 <- c(
    cos(omega / 2),
    sin(omega / 2) * c(
      sin(theta) * cos(phi),
      sin(theta) * sin(phi),
      cos(theta)
    )
  )
  q2list <- purrr::pmap(list(qts2$w, qts2$x, qts2$y, qts2$z), c)
  q2list <- purrr::map(q2list, ~ multiply_quaternions(q0, .x))
  qts2$w <- purrr::map_dbl(q2list, 1)
  qts2$x <- purrr::map_dbl(q2list, 2)
  qts2$y <- purrr::map_dbl(q2list, 3)
  qts2$z <- purrr::map_dbl(q2list, 4)
  DTW(
    qts1 = qts1,
    qts2 = qts2,
    resample = resample,
    disable_normalization = disable_normalization,
    distance_only = distance_only,
    step_pattern = step_pattern
  )
}

multiply_quaternions <- function(q1, q2) {
  c(
    q1[1] * q2[1] - sum(q1[-1] * q2[-1]),
    q1[1] * q2[2] + q1[2] * q2[1] - q1[3] * q2[4] - q1[4] * q2[3],
    q1[1] * q2[3] + q1[3] * q2[1] + q1[4] * q2[2] - q1[2] * q2[4],
    q1[1] * q2[4] + q1[4] * q2[1] + q1[2] * q2[3] - q1[3] * q2[2]
  )
}
