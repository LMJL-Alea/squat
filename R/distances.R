#' Distances between quaternions
#'
#' This is a collection of functions that assess the distance between two
#' quaternions in different spaces.
#'
#' @param q1 A \code{\link[onion]{quaternion}} object.
#' @param q2 A \code{\link[onion]{quaternion}} object.
#'
#' @return A positive scalar evaluating the distance between the two input
#'   quaternions.
#' @name distances
#'
#' @examples
#' q1 <- onion::rquat(1, rand = "norm")
#' q2 <- onion::rquat(1, rand = "norm")
#' geodesic(q1, q2)
NULL

#' @export
#' @rdname distances
geodesic <- function(q1, q2) {
  if (round(q1, 3) == round(q2, 3))
    return(0)
  2 * acos(Re(Conj(q1) * q2))
}
