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
#' @return A time series stored as a [tibble::tibble] with columns
#'   `time` and `distance` in which `distance` measures the angular distance
#'   between the quaternions of both input QTS at a given time point.
#'
#' @export
#' @examples
#' qts2distance(vespa64$igp[[1]], vespa64$igp[[2]])
qts2dts <- function(x, y) {
  if (!is_qts(x))
    cli::cli_abort("The input argument {.arg x} should be of class {.cls qts}.")
  if (!is_qts(y))
    cli::cli_abort("The input argument {.arg y} should be of class {.cls qts}.")
  if (!all(x$time == y$time))
    cli::cli_abort("The two input QTS should be evaluated on the same time grid.")
  qts2dts_impl(x, y)
}
