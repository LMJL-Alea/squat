#' Distance Matrix for Quaternion Time Series Samples
#'
#' @param t A list of time grids over which the QTS are measured.
#' @param q A list of QTS.
#' @param labels A character vector specifying labels for each QTS.
#'
#' @return A \code{\link[stats]{dist}} object storing the distance matrix
#'   between QTS in a sample via DTW.
#' @export
#'
#' @examples
#' # TO DO
distDTW <- function(t, q, labels = NULL) {
  n <- length(t)
  if (is.null(labels)) labels <- 1:n
  d <- numeric(n * (n - 1) / 2)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      d[n*(i-1) - i*(i-1)/2 + j-i] <- DTW(
        s1 = q[[i]], s2 = q[[j]],
        t1 = t[[i]], t2 = t[[j]],
        distance_only = TRUE
      )$normalizedDistance
    }
  }
  attributes(d) <- NULL
  attr(d, "Labels") <- labels
  attr(d, "Size") <- n
  attr(d, "call") <- match.call()
  class(d) <- "dist"
  d
}
