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
#' @name dtw
#'
#' @examples
#' # TO DO
NULL

#' @export
#' @rdname dtw
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

#' @export
#' @rdname dtw
DTWi <- function(qts1, qts2,
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

  lb <- rep(0, 6)
  ub <- rep(c(pi, 2*pi, 2*pi), times = 2)

  # Find suitable initial guess by global optim
  opt <- nloptr::directL(
    fn = cost,
    lower = lb,
    upper = ub,
    qts1 = qts1, qts2 = qts2,
    resample = FALSE,
    disable_normalization = TRUE,
    distance_only = TRUE,
    step_pattern = step_pattern
  )

  par <- opt$par
  for (i in 1:6) {
    if (par[i] < 0)
      par[i] <- .Machine$double.eps^0.5
  }

  if (par[1] > pi)
    par[1] <- pi - .Machine$double.eps^0.5

  if (par[2] > 2*pi)
    par[2] <- 2*pi - .Machine$double.eps^0.5

  if (par[3] > 2*pi)
    par[3] <- 2*pi - .Machine$double.eps^0.5

  if (par[4] > pi)
    par[4] <- pi - .Machine$double.eps^0.5

  if (par[5] > 2*pi)
    par[5] <- 2*pi - .Machine$double.eps^0.5

  if (par[6] > 2*pi)
    par[6] <- 2*pi - .Machine$double.eps^0.5

  # Run local optim
  opt <- stats::optim(
    par = par,
    fn = cost,
    lower = lb,
    upper = ub,
    method = "L-BFGS-B",
    qts1 = qts1, qts2 = qts2,
    resample = FALSE,
    disable_normalization = TRUE,
    distance_only = TRUE,
    step_pattern = step_pattern
  )

  res <- cost_impl(
    x = opt$par,
    qts1 = qts1, qts2 = qts2,
    resample = FALSE,
    disable_normalization = TRUE,
    distance_only = distance_only,
    step_pattern = step_pattern
  )
  res$par <- opt$par
  res
}

cost <- function(x,
                 qts1, qts2,
                 resample,
                 disable_normalization,
                 distance_only,
                 step_pattern) {
  cost_impl(
    x,
    qts1, qts2,
    resample,
    disable_normalization,
    distance_only,
    step_pattern
  )$distance
}

cost_impl <- function(x,
                      qts1, qts2,
                      resample,
                      disable_normalization,
                      distance_only,
                      step_pattern) {
  qts1 <- reorient_from_angles(qts1, x[1], x[2], x[3])
  qts2 <- reorient_from_angles(qts2, x[4], x[5], x[6])

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
    q1[1] * q2[2] + q1[2] * q2[1] + q1[3] * q2[4] - q1[4] * q2[3],
    q1[1] * q2[3] + q1[3] * q2[1] + q1[4] * q2[2] - q1[2] * q2[4],
    q1[1] * q2[4] + q1[4] * q2[1] + q1[2] * q2[3] - q1[3] * q2[2]
  )
}

reorient_from_angles <- function(qts, theta, phi, omega) {
  sgn <- (qts$w[1] > 0)
  q0 <- c(
    cos(omega / 2),
    sin(omega / 2) * c(
      sin(theta) * cos(phi),
      sin(theta) * sin(phi),
      cos(theta)
    )
  )

  qlist <- purrr::pmap(list(qts$w, qts$x, qts$y, qts$z), c)
  qlist <- purrr::map(qlist, ~ multiply_quaternions(q0, .x))
  qts$w <- purrr::map_dbl(qlist, 1)
  qts$x <- purrr::map_dbl(qlist, 2)
  qts$y <- purrr::map_dbl(qlist, 3)
  qts$z <- purrr::map_dbl(qlist, 4)

  if (xor(qts$w[1] > 0, sgn))
    qts[, 2:5] <- -qts[, 2:5]

  qts
}
