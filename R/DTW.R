#' Dynamic Time Warping for Quaternion Time Series
#'
#' This function evaluates the Dynamic Time Warping (DTW) distance between two
#' quaternion time series (QTS).
#'
#' If no evaluation grid is provided, the function assumes that the two input
#' QTS are evaluated on the same grid.
#'
#' @param q1 A \code{4 x p1} matrix representing the first QTS.
#' @param q2 A \code{4 x p2} matrix representing the second QTS.
#' @param t1 An optional numeric vector of size \code{p1} specifying the
#'   evaluation grid for the first QTS.
#' @param t2 An optional numeric vector of size \code{p2} specifying the
#'   evaluation grid for the second QTS.
#' @param distance_only A boolean specifyung whether to only compute distance
#'   (no backtrack, faster). Default is \code{FALSE}.
#'
#' @return An object of class \code{\link[dtw]{dtw}} storing the dynamic time
#'   warping results.
#' @export
#'
#' @examples
#' s1_raw <- onion::rquat(15, rand = "norm")
#' s1 <- s1_raw / Mod(s1_raw)
#' s2_raw <- onion::rquat(20, rand = "norm")
#' s2 <- s2_raw / Mod(s2_raw)
#' t1 <- seq(0, 1, length.out = 15)
#' t2 <- seq(0, 1, length.out = 20)
#' DTW(s1, s2, t1, t2)
DTW <- function(s1, s2, t1 = NULL, t2 = NULL, distance_only = FALSE) {
  s1 <- as.matrix(s1)
  s2 <- as.matrix(s2)
  if (!is.null(t1))
    s1 <- RegularizeGrid(t1, s1)
  if (!is.null(t2))
    s2 <- RegularizeGrid(t2, s2)
  M <- GetCostMatrix(s1, s2)
  dtw::dtw(M, distance.only = distance_only)
}
