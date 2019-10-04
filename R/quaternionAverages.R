#' Fréchet Means for Quaternions
#'
#' This is a collection of functions that provide access to Fréchet means of
#' quaternion samples according to various distances.
#'
#' @param x An object of class \code{\link[onion]{quaterion}} storing a sample
#'   of quaternions.
#'
#' @return An object of class \code{\link[onion]{quaterion}} and size 1 storing
#'   the sample mean of the input sample of quaternions.
#' @name averages
#'
#' @examples
#' qdata <- onion::rquat(15, rand = "norm")
#' mean(qdata)
NULL

#' @export
#' @rdname averages
mean.quaternion <- function(x) {
  M <- x %*% t(x)
  onion::as.quaternion(svd(M)$u[, 1, drop = FALSE])
}
