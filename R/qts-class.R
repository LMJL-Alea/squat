#' QTS Class
#'
#' A collection of functions that implements the QTS class. It currently
#' provides the `as_qts()` function for QTS coercion of [tibble::tibble]s and
#' the `is_qts()` function for checking if an object is a QTS.
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

derivative_qts <- function(qts) {
  qts <- derivative_qts_impl(qts)
  as_qts(qts[-1, ])
}

straighten_qts <- function(qts) {
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

log_qts <- function(x) {
  x <- log_qts_impl(x)
  as_qts(x)
}

exp_qts <- function(x) {
  x <- exp_qts_impl(x)
  as_qts(x)
}

reorient_qts <- function(x, disable_normalization = FALSE) {
  if (!disable_normalization) x <- normalize_qts(x)
  x <- reorient_qts_impl(x)
  as_qts(x)
}

normalize_qts <- function(x) {
  x <- normalize_qts_impl(x)
  as_qts(x)
}

#' QTS Centering and Standardization
#'
#' This function operates a centering of the QTS around the geometric mean of
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
#' centring(vespa64$igp[[1]])
centring <- function(x, standardize = FALSE, keep_summary_stats = FALSE) {
  if (!is_qts(x))
    cli::cli_abort("The input argument {.arg x} should be of class {.cls qts}.")
  out <- centring_qts_impl(x, standardize = standardize)
  out$qts <- as_qts(out$qts)
  if (keep_summary_stats) return(out)
  out$qts
}

resample_qts <- function(x,
                         tmin = NA, tmax = NA, nout = 0L,
                         disable_normalization = FALSE) {
  if (!disable_normalization)
    x <- normalize_qts(x)
  x <- resample_qts_impl(x, tmin, tmax, nout)
  as_qts(x)
}

smooth_qts <- function(x, alpha = 0.5) {
  x <- smooth_qts_impl(x, alpha)
  as_qts(x)
}

#' QTS Visualization
#'
#' @param x An object of class [qts].
#' @param change_points An integer vector specifying the indices of the change
#'   points to display if any. Defaults to `NULL`, in which case no change
#'   points are displayed.
#' @param ... Further arguments to be passed to methods.
#'
#' @return The [plot.qts()] method does not return anything while the
#'   [autoplot.qts()] method returns a [ggplot2::ggplot] object.
#'
#' @importFrom graphics plot
#' @export
#'
#' @examples
#' plot(vespa64$igp[[1]])
#' ggplot2::autoplot(vespa64$igp[[1]])
plot.qts <- function(x, change_points = NULL, ...) {
  print(autoplot(x, change_points = change_points, ...))
}

#' @importFrom ggplot2 autoplot .data
#' @export
#' @rdname plot.qts
autoplot.qts <- function(x, change_points = NULL, ...) {
  if (!is.null(change_points)) {
    if (!all(change_points %in% 1:nrow(x)))
      cli::cli_abort("The change point indices are out of bounds.")
    change_points <- x$time[change_points]
  }
  x <- tidyr::pivot_longer(x, cols = "w":"z")
  p <- ggplot2::ggplot(x, ggplot2::aes(
    x = .data$time,
    y = .data$value
  )) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(ggplot2::vars(.data$name), ncol = 1, scales = "free") +
    ggplot2::theme_linedraw() +
    ggplot2::labs(
      title = "Quaternion Time Series",
      x = "Time",
      y = ""
    )

  if (!is.null(change_points)) {
    p <- p +
      ggplot2::geom_point(
        data = subset(x, x$time %in% change_points),
        color = "red"
      )
  }

  p
}
