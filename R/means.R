#' Fréchet Means for Quaternions
#'
#' This is a collection of functions that provide access to Fréchet means of
#' quaternion samples according to various distances.
#'
#' @param x An object of class \code{\link[onion]{quaternion}} storing a sample
#'   of quaternions.
#'
#' @return An object of class \code{\link[onion]{quaternion}} and size 1 storing
#'   the sample mean of the input sample of quaternions.
#' @name averages
#'
#' @examples
#' qdata <- onion::rquat(15, rand = "norm")
#' mean(qdata)
NULL

#' @export
#' @rdname averages
mean_qts <- function(q, t = NULL) {
  grid_size <- 101
  # Regularise if not already done
  if (!is.null(t))
    q <- purrr::map2(t, q, RegularizeGrid, outSize = grid_size)
  1:grid_size %>%
    purrr::map(~ t(sapply(q, function(.y) .y[, .x]))) %>%
    sapply(function(.x) {
    .x %>%
      rotations::as.Q4() %>%
      mean(type = "geometric")
  })
}
