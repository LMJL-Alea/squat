#' Unit Quaternion Geodesic Distance
#'
#' This function computes the geodesic distance between two unit quaternions.
#'
#' @param x A length-4 numeric vector of unit norm representing the first
#'   quaternion.
#' @param y A length-4 numeric vector of unit norm representing the second
#'   quaternion.
#'
#' @return A positive scalar providing a measure of distance between the two
#'   input quaternions.
#'
#' @export
gdistance <- function(x, y) {
  x <- abs(sum(x * y))
  if (x > 1) x <- 1
  2 * acos(x)
}

qinv <- function(x) {
  x[2:4] <- -x[2:4]
  x
}

qprod <- function(x, y) {
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  x4 <- x[4]
  y1 <- y[1]
  y2 <- y[2]
  y3 <- y[3]
  y4 <- y[4]
  c(
    x1 * y1 - x2 * y2 - x3 * y3 - x4 * y4,
    x1 * y2 + x2 * y1 + x3 * y4 - x4 * y3,
    x1 * y3 + x3 * y1 - x2 * y4 + x4 * y2,
    x4 * y1 + x1 * y4 + x2 * y3 - x3 * y2
  )
}

#' Quaternion time series normalization
#'
#' @param x A 4-row numeric matrix.
#'
#' @return A 4-row numeric matrix.
#'
#' @export
#' @examples
#' qts <- matrix(rnorm(40), 4, 10)
#' qts <- apply(qts, 2, function(.x) .x / sqrt(sum(.x^2)))
#' normalize_qts(qts)
normalize_qts <- function(x, normalizer = NULL) {
  if (is.null(normalizer)) normalizer <- x[, 1]
  q_init <- qinv(normalizer)
  apply(x, 2, function(.x) qprod(q_init, .x))
}
