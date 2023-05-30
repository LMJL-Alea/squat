#' QTS K-Means Alignment Algorithm
#'
#' This function massages the input quaternion time series to feed them into the
#' k-means alignment algorithm for jointly clustering and aligning the input
#' QTS.
#'
#' @param x Either a numeric matrix of data, or an object that can be coerced to
#'   such a matrix (such as a numeric vector or a data frame with all numeric
#'   columns) or an object of class [qts_sample].
#' @param n_clusters An integer value specifying the number of clusters to be
#'   look for.
#' @param iter_max An integer value specifying the maximum number of iterations
#'   for terminating the k-mean algorithm. Defaults to `10L`.
#' @inheritParams stats::kmeans
#' @inheritParams fdacluster::fdakmeans
#'
#' @return An object of class [`stats::kmeans`] if the input `x` is NOT of class
#'   [`qts_sample`]. Otherwise, an object of class `kma_qts` which is
#'   effectively a list with three components:
#' - `qts_aligned`: An object of class [qts_sample] storing the sample of
#' aligned QTS;
#' - `qts_centers`: A list of objects of class [qts] representing the centers of
#' the clusters;
#' - `best_kma_result`: An object of class [fdacluster::caps] storing the
#' results of the best k-mean alignment result among all initialization that
#' were tried.
#'
#' @export
#' @examples
#' res_kma <- kmeans(vespa64$igp[1:10], n_clusters = 2)
kmeans <- function(x, n_clusters, ...) {
  UseMethod("kmeans")
}

#' @export
#' @rdname kmeans
kmeans.default <- function(x,
                           n_clusters = 1,
                           iter_max = 10,
                           nstart = 1,
                           algorithm =  c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"),
                           trace = FALSE,
                           ...) {
  stats::kmeans(
    x = x,
    centers = n_clusters,
    iter.max = iter_max,
    nstart = nstart,
    algorithm = algorithm,
    trace = trace
  )
}

#' @export
#' @rdname kmeans
kmeans.qts_sample <-function(x,
                             n_clusters = 1L,
                             seeds = NULL,
                             seeding_strategy = c("kmeans++", "exhaustive-kmeans++", "exhaustive", "hclust"),
                             warping_class = c("affine", "dilation", "none", "shift", "srsf"),
                             centroid_type = "mean",
                             metric = c("l2", "pearson"),
                             cluster_on_phase = FALSE,
                             ...) {
  if (!is_qts_sample(x))
    cli::cli_abort("The input argument {.arg x} should be of class {.cls qts_sample}.")

  q_list <- log(x)
  q_list <- purrr::map(q_list, \(.x) rbind(.x$x, .x$y, .x$z))
  t_list <- purrr::map(x, "time")

  # Prep data
  N <- length(q_list)
  L <- dim(q_list[[1]])[1]
  P <- dim(q_list[[1]])[2]

  if (is.null(t_list))
    grid <- 0:(P-1)
  else
    grid <- matrix(nrow = N, ncol = P)

  values <- array(dim = c(N, L, P))
  for (n in 1:N) {
    values[n, , ] <- q_list[[n]]
    if (!is.null(t_list)) {
      grid[n, ] <- t_list[[n]]
    }
  }

  out <- fdacluster::fdakmeans(
    x = grid,
    y = values,
    n_clusters = n_clusters,
    seeds = seeds,
    seeding_strategy = seeding_strategy,
    warping_class = warping_class,
    centroid_type = centroid_type,
    metric = metric,
    cluster_on_phase = cluster_on_phase
  )

  res <- list(
    qts_aligned = purrr::imap(out$memberships, \(.label, .id) {
      exp(as_qts(tibble::tibble(
        time = out$grids[.label, ],
        w    = 0,
        x    = out$aligned_curves[.id, 1, ],
        y    = out$aligned_curves[.id, 2, ],
        z    = out$aligned_curves[.id, 3, ]
      )))
    }),
    qts_centers = purrr::map(1:out$n_clusters, \(.id) {
      exp(as_qts(tibble::tibble(
        time = out$grids[.id, ],
        w    = 0,
        x    = out$center_curves[.id, 1, ],
        y    = out$center_curves[.id, 2, ],
        z    = out$center_curves[.id, 3, ]
      )))
    }),
    best_kma_result = out
  )

  class(res) <- "kma_qts"
  res
}

#' Plot for `kma_qts` objects
#'
#' This function creates a visualization of the results of the k-means alignment
#' algorithm applied on a sample of QTS and returns the corresponding
#' [ggplot2::ggplot] object which enable further customization of the plot.
#'
#' @param object An object of class `kma_qts` as produced by the
#'   [kmeans.qts_sample()] method.
#' @param ... Further arguments to be passed to other methods.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @importFrom ggplot2 autoplot .data
#' @export
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' res_kma <- kmeans(vespa64$igp[1:10], k = 2, nstart = 1)
#' ggplot2::autoplot(res_kma)
autoplot.kma_qts <- function(object, ...) {
  data <- as_qts_sample(c(object$qts_centers, object$qts_aligned))
  n <- length(object$qts_aligned)
  k <- length(object$qts_centers)
  memb <- c(1:k, object$best_kma_result$memberships)
  high <- c(rep(TRUE, k), rep(FALSE, n))
  autoplot(data, memberships = memb, highlighted = high) +
    ggplot2::labs(
      title = "K-Means Alignment Clustering Results",
      subtitle = cli::pluralize("Using {k} cluster{?s}")
    )
}

#' Plot for `kma_qts` objects
#'
#' This function creates a visualization of the results of the k-means alignment
#' algorithm applied on a sample of QTS **without** returning the plot data as
#' an object.
#'
#' @param x An object of class `kma_qts` as produced by the [kmeans()] function.
#' @inheritParams autoplot.kma_qts
#'
#' @return NULL
#'
#' @importFrom graphics plot
#' @export
#' @examples
#' res_kma <- kmeans(vespa64$igp[1:10], k = 2, nstart = 1)
#' plot(res_kma)
plot.kma_qts <- function(x, ...) {
  print(autoplot(x, ...))
}
