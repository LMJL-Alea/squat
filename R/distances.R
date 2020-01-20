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
