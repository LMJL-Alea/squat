#' QTS Transformation To Distance Time Series
#'
#' This function computes a real-valued time series reporting the pointwise
#' geodesic distance between the two input QTS at each time point.
#'
#' The function currently expects that the two input QTS are evaluated on the
#' same time grid.
#'
#' @param x An object of class [qts].
#' @param y An object of class [qts].
#'
#' @return A time series stored as a [tibble::tibble] with columns `time` and
#'   `distance` in which `distance` measures the angular distance between the
#'   quaternions of both input QTS at a given time point.
#'
#' @export
#' @examples
#' qts2dts(vespa64$igp[[1]], vespa64$igp[[2]])
qts2dts <- function(x, y) {
  if (!is_qts(x)) {
    cli::cli_abort("The input argument {.arg x} should be of class {.cls qts}.")
  }
  if (!is_qts(y)) {
    cli::cli_abort("The input argument {.arg y} should be of class {.cls qts}.")
  }
  if (!all(x$time == y$time)) {
    cli::cli_abort(
      "The two input QTS should be evaluated on the same time grid."
    )
  }
  qts2dts_impl(x, y)
}

#' QTS Transformation To Norm Time Series
#'
#' This function computes a univariate time series representing the norm of the
#' quaternions.
#'
#' @param x An object of class [qts].
#' @param disable_normalization A boolean specifying whether quaternion
#'   normalization should be disabled. Defaults to `FALSE`.
#'
#' @return A time series stored as a [tibble::tibble] with columns `time` and
#'   `norm` in which `norm` measures the angular distance between the current
#'   quaternion and the identity.
#'
#' @export
#' @examples
#' qts2nts(vespa64$igp[[1]])
qts2nts <- function(x, disable_normalization = FALSE) {
  if (!is_qts(x)) {
    cli::cli_abort("The input argument {.arg x} should be of class {.cls qts}.")
  }
  qts2nts_impl(x, disable_normalization = disable_normalization)
}

#' QTS Transformation To Angle Time Series
#'
#' This function computes a univariate time series representing the angle
#' between the first and other attitudes.
#'
#' @param x An object of class [qts].
#' @param disable_normalization A boolean specifying whether quaternion
#'   normalization should be disabled. Defaults to `FALSE`.
#'
#' @return A time series stored as a [tibble::tibble] with columns `time` and
#'   `angle` in which `angle` measures the angle between the current rotation
#'   and the first one.
#'
#' @export
#' @examples
#' qts2ats(vespa64$igp[[1]])
qts2ats <- function(x, disable_normalization = FALSE) {
  if (!is_qts(x)) {
    cli::cli_abort("The input argument {.arg x} should be of class {.cls qts}.")
  }
  out <- qts2ats_impl(x, disable_normalization = disable_normalization)
  out$angle[out$angle < .Machine$double.eps] <- 0
  out
}

#' QTS Transformation to Angular Velocity Vector Time Series
#'
#' This function projects a quaternion time series into the space of angular
#' velocity vectors.
#'
#' @param x An object of class [qts].
#' @inheritParams stats::smooth.spline
#'
#' @return A time series stored as a [tibble::tibble] with columns `time`, `x`,
#'   `y` and `z` containing the angular velocity vector at each time point.
#'
#' @export
#' @examples
#' qts2avvts(vespa64$igp[[1]])
qts2avvts <- function(x, spar = 0) {
  if (!is_qts(x)) {
    cli::cli_abort("The input argument {.arg x} should be of class {.cls qts}.")
  }
  tmp <- log(x)
  mod_x <- stats::smooth.spline(tmp$time, tmp$x, spar = spar)
  mod_y <- stats::smooth.spline(tmp$time, tmp$y, spar = spar)
  mod_z <- stats::smooth.spline(tmp$time, tmp$z, spar = spar)
  tibble::tibble(
    time = tmp$time,
    x = 2 * stats::predict(mod_x, tmp$time, deriv = 1)$y,
    y = 2 * stats::predict(mod_y, tmp$time, deriv = 1)$y,
    z = 2 * stats::predict(mod_z, tmp$time, deriv = 1)$y
  )
}

#' QTS Transformation to Angular Velocity Magnitude Time Series
#'
#' This function projects a quaternion time series into the space of angular
#' velocity magnitudes.
#'
#' @param x An object of class [qts].
#' @inheritParams stats::smooth.spline
#'
#' @return A time series stored as a [tibble::tibble] with columns `time` and
#'   `magnitude` containing the angular velocity magnitude at each time point.
#'
#' @export
#' @examples
#' qts2avmts(vespa64$igp[[1]])
qts2avmts <- function(x, spar = 0) {
  out <- qts2avvts(x, spar = spar)
  out$magnitude <- sqrt(out$x^2 + out$y^2 + out$z^2)
  out$x <- out$y <- out$z <- NULL
  out
}

#' QTS Transformation to Smoothed Quaternion Time Series
#'
#' This function smooths a given QTS using a b-spline functional representation.
#'
#' @param x An object of class [qts].
#' @inheritParams stats::smooth.spline
#'
#' @return An object of class `qts` storing the smoothed QTS.
#'
#' @export
#' @examples
#' qts2sqts(vespa64$igp[[1]], spar = 0.3)
qts2sqts <- function(x, spar = 0) {
  if (!is_qts(x)) {
    cli::cli_abort("The input argument {.arg x} should be of class {.cls qts}.")
  }
  tmp <- log(x)
  exp(as_qts(tibble::tibble(
    time = unique(tmp$time),
    w = 0,
    x = stats::smooth.spline(tmp$time, tmp$x, spar = spar)$y,
    y = stats::smooth.spline(tmp$time, tmp$x, spar = spar)$y,
    z = stats::smooth.spline(tmp$time, tmp$x, spar = spar)$y
  )))
}

#' QTS Transformation to Angular Acceleration Vector Time Series
#'
#' This function projects a quaternion time series into the space of angular
#' acceleration vectors.
#'
#' @param x An object of class [qts].
#' @inheritParams stats::smooth.spline
#'
#' @return A time series stored as a [tibble::tibble] with columns `time`, `x`,
#'   `y` and `z` containing the angular acceleration vector at each time point.
#'
#' @export
#' @examples
#' qts2aavts(vespa64$igp[[1]])
qts2aavts <- function(x, spar = 0) {
  if (!is_qts(x)) {
    cli::cli_abort("The input argument {.arg x} should be of class {.cls qts}.")
  }
  tmp <- log(x)
  mod_x <- stats::smooth.spline(tmp$time, tmp$x, spar = spar)
  mod_y <- stats::smooth.spline(tmp$time, tmp$y, spar = spar)
  mod_z <- stats::smooth.spline(tmp$time, tmp$z, spar = spar)
  tibble::tibble(
    time = tmp$time,
    x = 2 * stats::predict(mod_x, tmp$time, deriv = 2)$y,
    y = 2 * stats::predict(mod_y, tmp$time, deriv = 2)$y,
    z = 2 * stats::predict(mod_z, tmp$time, deriv = 2)$y
  )
}

#' QTS Transformation to Angular Acceleration Magnitude Time Series
#'
#' This function projects a quaternion time series into the space of angular
#' acceleration magnitudes.
#'
#' @param x An object of class [qts].
#' @inheritParams stats::smooth.spline
#'
#' @return A time series stored as a [tibble::tibble] with columns `time` and
#'   `mag` containing the angular acceleration magnitude at each time point.
#'
#' @export
#' @examples
#' qts2aamts(vespa64$igp[[1]])
qts2aamts <- function(x, spar = 0) {
  out <- qts2aavts(x, spar = spar)
  out$magnitude <- sqrt(out$x^2 + out$y^2 + out$z^2)
  out$x <- out$y <- out$z <- NULL
  out
}

#' QTS Transformation to Angle-Axis Time Series
#'
#' This function converts a quaternion time series into its angle-axis
#' representation.
#'
#' @param x An object of class [qts].
#'
#' @return A time series stored as a [tibble::tibble] with columns `time`,
#'   `angle`, `ux`, `uy` and `uz` containing the angle-axis representation of
#'   the input quaternions.
#'
#' @export
#' @examples
#' qts2aats(vespa64$igp[[1]])
qts2aats <- function(x) {
  if (!is_qts(x)) {
    cli::cli_abort("The input argument {.arg x} should be of class {.cls qts}.")
  }
  out <- qts2aats_impl(x)
  names(out) <- c("time", "angle", "ux", "uy", "uz")
  out
}

#' QTS Transformation to Roll-Pitch-Yaw Time Series
#'
#' This function converts a quaternion time series into its roll-pitch-yaw
#' angles representation.
#'
#' @param x An object of class [qts].
#'
#' @return A time series stored as a [tibble::tibble] with columns `time`,
#'   `roll`, `pitch` and `yaw` containing the roll, pitch and yaw angles
#'   representation of the input quaternions.
#'
#' @export
#' @examples
#' qts2rpyts(vespa64$igp[[1]])
qts2rpyts <- function(x) {
  if (!is_qts(x)) {
    cli::cli_abort("The input argument {.arg x} should be of class {.cls qts}.")
  }
  out <- qts2rpyts_impl(x)
  names(out) <- c("time", "roll", "pitch", "yaw")
  out
}

rpyts2qts <- function(x) {
  x |>
    rpyts2qts_impl() |>
    as_qts()
}
