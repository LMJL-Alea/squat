#' Distance Matrix for Quaternion Time Series Samples
#'
#' @param q A list of QTS.
#' @param t An optional list of the same size as \code{q} containing evaluation
#'   grids for each QTS.
#' @param labels A character vector specifying labels for each QTS.
#' @param step_pattern The choice for the step pattern of the warping path (see
#'  \code{\link[dtw]{stepPattern}} for more information)
#'
#' @return A \code{\link[stats]{dist}} object storing the distance matrix
#'   between QTS in a sample via DTW.
#' @export
#'
#' @examples
#' # TO DO

distDTW <- function (q, t = NULL, labels = NULL, step_pattern = dtw::symmetric2)
{
  n <- length(q)
  if (is.null(labels))
    labels <- 1:n
  d <- numeric(n * (n - 1)/2)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (is.null(t))
        d[n * (i - 1) - i * (i - 1)/2 + j - i] <- DTW(s1 = q[[i]],
                                                      s2 = q[[j]], distance_only = TRUE, step_pattern = step_pattern)$distance
      else d[n * (i - 1) - i * (i - 1)/2 + j - i] <- DTW(s1 = q[[i]],
                                                         s2 = q[[j]], t1 = t[[i]], t2 = t[[j]], distance_only = TRUE, step_pattern = step_pattern)$distance
    }
  }
  attributes(d) <- NULL
  attr(d, "Labels") <- labels
  attr(d, "Size") <- n
  attr(d, "call") <- match.call()
  class(d) <- "dist"
  d
}
