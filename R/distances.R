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
  x <- abs(1 - sum((x - y)^2) / 2)
  if (x > 1) x <- 1
  2 * acos(x)
}

#' Quaternion Inverse
#'
#' @param q A length-4 vector storing a quaternion.
#'
#' @return A length-4 vector storing the inverse of the input quaternion.
#' @export
#'
#' @examples
#' q <- rnorm(4)
#' q <- q / sqrt(sum(q^2))
#' q
#' qinv(q)
qinv <- function(q) {
  q[2:4] <- -q[2:4]
  q
}

#' Quaternion Product
#'
#' @param q1 A length-4 vector storing a quaternion.
#' @param q2 A length-4 vector storing a quaternion.
#'
#' @return A length-4 vector storing the product of the two input quaternions.
#' @export
#'
#' @examples
#' q1 <- rnorm(4)
#' q1 <- q1 / sqrt(sum(q1^2))
#' q2 <- rnorm(4)
#' q2 <- q2 / sqrt(sum(q2^2))
#' qprod(q1, q2)
qprod <- function(q1, q2) {
  q1 <- as.numeric(q1)
  q2 <- as.numeric(q2)
  x1 <- q1[1]
  x2 <- q1[2]
  x3 <- q1[3]
  x4 <- q1[4]
  y1 <- q2[1]
  y2 <- q2[2]
  y3 <- q2[3]
  y4 <- q2[4]
  c(
    w = x1 * y1 - x2 * y2 - x3 * y3 - x4 * y4,
    x = x1 * y2 + x2 * y1 + x3 * y4 - x4 * y3,
    y = x1 * y3 + x3 * y1 - x2 * y4 + x4 * y2,
    z = x4 * y1 + x1 * y4 + x2 * y3 - x3 * y2
  )
}

#' Quaternion time series normalization
#'
#' @param x A 4-row numeric matrix.
#' @param normalizer A length-4 vector storing a quaternion to use as origin.
#'   Defaults to \code{NULL}, which takes the first quaternion in QTS \code{x}
#'   as origin.
#'
#' @return A 4-row numeric matrix storing the normalized QTS.
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
