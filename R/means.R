#' Fréchet Means for Quaternions
#'
#' This is a collection of functions that provide access to Fréchet means of
#' quaternion samples according to various distances.
#'
#' @param q A list of QTS.
#' @param t An optional list of the same size as \code{q} containing evaluation
#'   grids for each QTS.
#'
#' @return An object of class \code{\link[onion]{quaternion}} and size 1 storing
#'   the sample mean of the input sample of quaternions.
#'
#' @export
#' @examples
#' qdata <- onion::rquat(15)
#' mean(qdata)
mean_qts <- function(q, t = NULL) {
  grid_size <- 101
  # Regularise if not already done
  if (!is.null(t))
    q <- purrr::map2(t, q, regularize_grid, out_size = grid_size)
  1:grid_size %>%
    purrr::map(~ t(sapply(q, function(.y) .y[, .x]))) %>%
    sapply(GetGeodesicMean)
}
