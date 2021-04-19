#' K-Means for Quaternion Time Series
#'
#' @param q A list of quaternion time series. Each QTS is a 4-row matrix.
#' @param t An optional list of evaluation grids for the QTS. By default, if you
#'   do not supply it, the function assumes that all QTS are evaluated on the
#'   same grid.
#' @param k The number of clusters (default: 2).
#' @param iter_max The maximum number of iterations before stopping the k-means
#'   algorithm if it failed converging (default: 20).
#' @param nstart The number of random restarts for the initial centers (default:
#'   1000).
#'
#' @return A list with two components: \code{memberships} is an integer vector
#'   specifying the membership of each data point and \code{Var} is the sum of
#'   within-sum-of-squares.
#' @export
#'
#' @examples
kmeans_qts <-function(q, t = NULL, k = 2, iter_max = 20, nstart = 1000) {
  # Assumes qsr is a list of matrices 4x101
  # and tsr is a vector of size 101
  n <- length(q)
  B <- choose(n, k)

  if (nstart > B)
    init <- utils::combn(n, k, simplify = FALSE)
  else
    init <- replicate(nstart, sample.int(n, k), simplify = FALSE)

  if (requireNamespace("furrr", quietly = TRUE)) {
    future::plan(future::multiprocess)
    solutions <- init %>%
      furrr::future_map(
        .f = kmeans_qts_single,
        q = q,
        t = t,
        iter_max = iter_max,
        .progress = TRUE
      )
  } else {
    solutions <- init %>%
      purrr::map(
        .f = kmeans_qts_single,
        q = q,
        t = t,
        iter_max = iter_max
      )
  }


  wss_vector <- purrr::map_dbl(solutions, "Var")
  solutions[[which.min(wss_vector)]]
}

kmeans_qts_single <- function(init, q, t = NULL, iter_max = 20) {
  # Step 0: initialization
  centers <- q[init]
  iter <- 0
  old_wss <- 0
  wss <- .Machine$double.xmax

  while (iter < iter_max & (wss < old_wss | iter == 0)) {
    old_wss <- wss

    # Step 1: Calculate distance to centers
    dist_to_centers <- q %>%
      purrr::map(function(.x) {
        centers %>%
          purrr::map_dbl(function(.y) {
            DTW(
              s1 = .x, s2 = .y,
              t1 = t, t2 = t,
              distance_only = TRUE
            )$normalizedDistance^2
          })
      })

    # Step 2: Update memberships
    memberships <- dist_to_centers %>%
      purrr::map_int(which.min)
    dist_to_centers <- dist_to_centers %>%
      purrr::map_dbl(min)

    # Step 3: Update centers
    centers <- centers %>%
      purrr::imap(~ mean_qts(q[memberships == .y]))

    wss <- sum(dist_to_centers)
    iter <- iter + 1
  }

  list(cluster = memberships, Var = mean(dist_to_centers))
}
