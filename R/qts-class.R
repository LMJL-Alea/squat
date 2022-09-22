#' QTS Class
#'
#' A collection of functions that implements the QTS class. It currently
#' provides the `as_qts()` function for QTS coercion of tibbles and the
#' `is_qts()` function for checking if an object is a QTS.
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
#' @param x A [tibble::tibble] with columns `time`, `w`, `x`, `y` and `z`.
#' @param digits An integer value specifying the number of digits to keep for
#'   printing. Defaults to `5L`.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An object of class [qts].
#' @name qts
#'
#' @examples
#' qts1 <- vespa64$igp[[1]]
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

#' @export
#' @rdname qts
format.qts <- function(x, digits = 5, ...) {
  x$w <- format_quaternion_component(x$w, digits = digits)
  x$x <- format_quaternion_component(x$x, digits = digits)
  x$y <- format_quaternion_component(x$y, digits = digits)
  x$z <- format_quaternion_component(x$z, digits = digits)
  NextMethod()
}

format_quaternion_component <- function(x, digits = 5) {
  tibble::num(
    x = round(x, digits = digits),
    digits = digits,
    notation = "dec"
  )
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
#' derivative_qts(vespa64$igp[[1]])
derivative_qts <- function(qts) {
  if (!is_qts(qts))
    cli::cli_abort("The input object should be of class {.cls qts}.")
  qts <- derivative_qts_impl(qts)
  as_qts(qts[-1, ])
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
#' straighten_qts(vespa64$igp[[1]])
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

#' QTS Logarithm
#'
#' This function computes the logarithm of a quaternion time series as the time
#' series of the quaternion logarithms.
#'
#' @param x An object of class \code{\link{qts}}.
#'
#' @return An object of class \code{\link{qts}} in which quaternions have been
#'   replaced by their logarithm.
#' @export
#'
#' @examples
#' log_qts(vespa64$igp[[1]])
log_qts <- function(x) {
  if (!is_qts(x))
    cli::cli_abort("The input argument {.arg x} should be of class {.cls qts}.")
  x <- log_qts_impl(x)
  as_qts(x)
}

#' QTS Exponential
#'
#' This function computes the exponential of a quaternion time series as the
#' time series of the quaternion exponentials.
#'
#' @param x An object of class \code{\link{qts}}.
#'
#' @return An object of class \code{\link{qts}} in which quaternions have been
#'   replaced by their exponential.
#' @export
#'
#' @examples
#' x <- log_qts(vespa64$igp[[1]])
#' exp_qts(x)
exp_qts <- function(x) {
  if (!is_qts(x))
    cli::cli_abort("The input argument {.arg x} should be of class {.cls qts}.")
  x <- exp_qts_impl(x)
  as_qts(x)
}

#' QTS Reorientation
#'
#' This function reorients the quaternions in a QTS for representing attitude
#' with respect to the orientation of the sensor at the first time point.
#'
#' @param x An object of class \code{\link{qts}}.
#' @param disable_normalization A boolean specifying whether quaternion
#'   normalization should be disabled. Defaults to `FALSE`.
#'
#' @return An object of class \code{\link{qts}} in which quaternions measure
#'   attitude with respect to the orientation of the sensor at the first time
#'   point.
#'
#' @export
#' @examples
#' reorient_qts(vespa64$igp[[1]])
reorient_qts <- function(x, disable_normalization = FALSE) {
  if (!is_qts(x))
    cli::cli_abort("The input argument {.arg x} should be of class {.cls qts}.")
  if (!disable_normalization) x <- normalize_qts(x)
  x <- reorient_qts_impl(x)
  as_qts(x)
}

#' QTS Normalization
#'
#' This function ensures that all quaternions in the time series are unit
#' quaternions.
#'
#' @param x An object of class \code{\link{qts}}.
#'
#' @return An object of class \code{\link{qts}} in which quaternions are unit
#'   quaternions.
#'
#' @export
#' @examples
#' normalize_qts(vespa64$igp[[1]])
normalize_qts <- function(x) {
  if (!is_qts(x))
    cli::cli_abort("The input argument {.arg x} should be of class {.cls qts}.")
  x <- normalize_qts_impl(x)
  as_qts(x)
}

#' QTS Centering and Standardization
#'
#' This function operates a centring of the QTS around the geometric mean of
#' its quaternions. This is effectively achieved by left-multiplying each
#' quaternion by the inverse of their geometric mean.
#'
#' @param x An object of class \code{\link{qts}}.
#' @param standardize A boolean specifying whether to standardize the QTS in
#'   addition to centering it. Defaults to `FALSE`.
#' @param keep_summary_stats A boolean specifying whether the mean and standard
#'   deviation used for standardizing the data should be stored in the output
#'   object. Defaults to `FALSE` in which case only the centered
#'   \code{\link{qts}} is returned.
#'
#' @return If `keep_summary_stats = FALSE`, an object of class \code{\link{qts}}
#'   in which quaternions have been centered (and possibly standardized) around
#'   their geometric mean. If `keep_summary_stats = TRUE`, a list with three
#'   components:
#'   - `qts`: an object of class \code{\link{qts}} in which quaternions have
#'   been centered (and possibly standardized) around their geometric mean;
#' - `mean`: a numeric vector with the quaternion Fréchet mean;
#' - `sd`: a numeric value with the quaternion Fréchet standard deviation.
#'
#' @export
#' @examples
#' centring_qts(vespa64$igp[[1]])
centring_qts <- function(x, standardize = FALSE, keep_summary_stats = FALSE) {
  if (!is_qts(x))
    cli::cli_abort("The input argument {.arg x} should be of class {.cls qts}.")
  out <- centring_qts_impl(x, standardize = standardize)
  out$qts <- as_qts(out$qts)
  if (keep_summary_stats) return(out)
  out$qts
}

#' QTS Resampling
#'
#' This function performs uniform resampling using SLERP.
#'
#' @param x An object of class [qts].
#' @param tmin A numeric value specifying the lower bound of the time interval
#'   over which uniform resampling should take place. It must satisfy `tmin >=
#'   min(qts$time)`. Defaults to `NA` in which case it is set to
#'   `min(qts$time)`.
#' @param tmax A numeric value specifying the upper bound of the time interval
#'   over which uniform resampling should take place. It must satisfy `tmax <=
#'   max(qts$time)`. Defaults to `NA` in which case it is set to
#'   `max(qts$time)`.
#' @param nout An integer specifying the size of the uniform grid for time
#'   resampling. Defaults to `0L` in which case it uses the same grid size as
#'   the input QTS.
#' @param disable_normalization A boolean specifying whether quaternion
#'   normalization should be disabled. Defaults to `FALSE` in which case the
#'   function makes sure that quaternions are normalized prior to performing
#'   SLERP interpolation.
#'
#' @return An object of class [qts] in which quaternions are uniformly sampled
#'   in the range `[tmin, tmax]`.
#'
#' @export
#' @examples
#' resample_qts(vespa64$igp[[1]])
resample_qts <- function(x,
                         tmin = NA, tmax = NA, nout = 0L,
                         disable_normalization = FALSE) {
  if (!is_qts(x))
    cli::cli_abort("The input argument {.arg x} should be of class {.cls qts}.")
  if (!disable_normalization)
    x <- normalize_qts(x)
  x <- resample_qts_impl(x, tmin, tmax, nout)
  as_qts(x)
}

#' QTS Smoothing via SLERP Interpolation
#'
#' This function performs a smoothing of a QTS by SLERP interpolation.
#'
#' @param x An object of class [qts].
#' @param alpha A numeric value in `[0,1]` specifying the amount of smoothing.
#'   The closer to one, the smoother the resulting QTS. Defaults to `0.5`.
#'
#' @return An object of class [qts] which is a smooth version of the input QTS.
#'
#' @export
#' @examples
#' smooth_qts(vespa64$igp[[1]])
smooth_qts <- function(x, alpha = 0.5) {
  if (!is_qts(x))
    cli::cli_abort("The input argument {.arg x} should be of class {.cls qts}.")
  x <- smooth_qts_impl(x, alpha)
  as_qts(x)
}
