#' Distance Matrix for Quaternion Time Series Samples
#'
#' @param q A list of QTS.
#' @param t An optional list of the same size as \code{q} containing evaluation
#'   grids for each QTS.
#' @param labels A character vector specifying labels for each QTS.
#' @inheritParams DTW
#' @param normalize A boolean specifying whether to compute normalized distance
#'   between QTS.
#'   Please note that not all step patterns are normalizable
#'   (default: \code{FALSE}).
#' @return A \code{\link[stats]{dist}} object storing the distance matrix
#'   between QTS in a sample via DTW.
#' @export
#'
#' @examples
#' # TO DO

distDTW <- function (q, t = NULL, labels = NULL, step_pattern = dtw::symmetric2, normalize = FALSE) {
  if (normalize && is.na(attr(step_pattern, "norm"))) stop("The provided step pattern is not normalizable.")
  n <- length(q)
  if (is.null(labels))
    labels <- 1:n
  d <- numeric(n * (n - 1)/2)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (is.null(t))
        if(normalize) {
          d[n * (i - 1) - i * (i - 1)/2 + j - i] <- DTW(
            s1 = q[[i]],
            s2 = q[[j]],
            distance_only = TRUE,
            step_pattern = step_pattern
          )$normalizedDistance
        }else{
          d[n * (i - 1) - i * (i - 1) / 2 + j - i] <- DTW(
          s1 = q[[i]],
          s2 = q[[j]],
          distance_only = TRUE,
          step_pattern = step_pattern
          )$distance
        }

      else{
        if(normalize) {
          d[n * (i - 1) - i * (i - 1)/2 + j - i] <- DTW(
            s1 = q[[i]],
            s2 = q[[j]],
            t1 = t[[i]],
            t2 = t[[j]],
            distance_only = TRUE,
            step_pattern = step_pattern
          )$normalizedDistance
        }else{
          d[n * (i - 1) - i * (i - 1)/2 + j - i] <- DTW(
            s1 = q[[i]],
            s2 = q[[j]],
            t1 = t[[i]],
            t2 = t[[j]],
            distance_only = TRUE,
            step_pattern = step_pattern
          )$distance
        }

      }
    }
  }
  attributes(d) <- NULL
  attr(d, "Labels") <- labels
  attr(d, "Size") <- n
  attr(d, "call") <- match.call()
  class(d) <- "dist"
  d
}
