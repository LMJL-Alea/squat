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
#' @param x A \code{\link[tibble]{tibble}} with columns `time`, `w`, `x`, `y`
#'   and `z`.
#'
#' @return An object of class \code{\link{qts}}.
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
  x$w <- tibble::num(x$w, notation = "dec")
  x$x <- tibble::num(x$x, notation = "dec")
  x$y <- tibble::num(x$y, notation = "dec")
  x$z <- tibble::num(x$z, notation = "dec")
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
