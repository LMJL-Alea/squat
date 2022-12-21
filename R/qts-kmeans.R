#' QTS K-Means Alignment Algorithm
#'
#' This function massages the input quaternion time series to feed them into the
#' k-means alignment algorithm for jointly clustering and aligning the input
#' QTS.
#'
#' @param x Either a numeric matrix of data, or an object that can be coerced to
#'   such a matrix (such as a numeric vector or a data frame with all numeric
#'   columns) or an object of class [qts_sample].
#' @param k An integer value specifying the number of clusters to be look for.
#' @param iter_max An integer value specifying the maximum number of iterations
#'   for terminating the k-mean algorithm. Defaults to `10L`.
#' @param nstart An integer value specifying the number of random restarts of
#'   the algorithm. The higher `nstart`, the more robust the result. Defaults to
#'   `1L`.
#' @inheritParams stats::kmeans
#' @param centroid A string specifying which type of centroid should be used
#'   when applying kmeans on a QTS sample. Choices are `mean` and `medoid`.
#'   Defaults to `mean`.
#' @param dissimilarity A string specifying which type of dissimilarity should
#'   be used when applying kmeans on a QTS sample. Choices are `l2` and
#'   `pearson`. Defaults to `l2`.
#' @param warping A string specifying which class of warping functions should be
#'   used when applying kmeans on a QTS sample. Choices are `none`, `shift`,
#'   `dilation` and `affine`. Defaults to `affine`.
#'
#' @return An object of class [`stats::kmeans`] if the input `x` is NOT of class
#'   [`qts_sample`]. Otherwise, an object of class `kma_qts` which is
#'   effectively a list with three components:
#' - `qts_aligned`: An object of class [qts_sample] storing the sample of
#' aligned QTS;
#' - `qts_centers`: A list of objects of class [qts] representing the centers of
#' the clusters;
#' - `best_kma_result`: An object of class [fdacluster::kma] storing the results
#' of the best k-mean alignment result among all initialization that were tried.
#'
#' @export
#' @examples
#' res_kma <- kmeans(vespa64$igp, k = 2)
kmeans <- function(x, k, iter_max = 10, nstart = 1, ...) {
  UseMethod("kmeans")
}

#' @export
#' @rdname kmeans
kmeans.default <- function(x,
                           k,
                           iter_max = 10,
                           nstart = 1,
                           algorithm =  c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"),
                           trace = FALSE,
                           ...) {
  stats::kmeans(
    x = x,
    centers = k,
    iter.max = iter_max,
    nstart = nstart,
    algorithm = algorithm,
    trace = trace
  )
}

#' @export
#' @rdname kmeans
kmeans.qts_sample <-function(x,
                             k = 1,
                             iter_max = 10,
                             nstart = 1,
                             centroid = "mean",
                             dissimilarity = "l2",
                             warping = "affine",
                             ...) {
  if (!is_qts_sample(x))
    cli::cli_abort("The input argument {.arg x} should be of class {.cls qts_sample}.")

  q_list <- purrr::map(x, log_qts)
  q_list <- purrr::map(q_list, ~ rbind(.x$x, .x$y, .x$z))
  t_list <- purrr::map(x, "time")

  # Prep data
  n <- length(q_list)
  d <- dim(q_list[[1]])[1]
  p <- dim(q_list[[1]])[2]

  if (is.null(t_list))
    grid <- 0:(p-1)
  else
    grid <- matrix(nrow = n, ncol = p)

  values <- array(dim = c(n, d, p))
  for (i in 1:n) {
    values[i, , ] <- q_list[[i]]
    if (!is.null(t_list)) {
      grid[i, ] <- t_list[[i]]
    }
  }

  B <- choose(n, k)

  if (nstart > B)
    init <- utils::combn(n, k, simplify = FALSE)
  else
    init <- replicate(nstart, sample.int(n, k), simplify = FALSE)

  .kma_with_multiple_restarts <- function(init_values) {
    pb <- progressr::progressor(along = init_values)
    furrr::future_map(init_values, ~ {
      pb()
      fdacluster::kma(
        x = grid,
        y = values,
        seeds = .x,
        n_clust = k,
        center_method = centroid,
        warping_method = warping,
        dissimilarity_method = dissimilarity,
        use_verbose = FALSE,
        warping_options = c(0.1, 0.1),
        use_fence = FALSE
      )
    }, .options = furrr::furrr_options(seed = TRUE))
  }

  solutions <- .kma_with_multiple_restarts(init)
  wss_vector <- purrr::map_dbl(solutions, ~ sum(.x$final_dissimilarity))
  opt <- solutions[[which.min(wss_vector)]]

  centers <- purrr::map(1:opt$n_clust_final, ~ {
    exp_qts(as_qts(tibble::tibble(
      time = opt$x_centers_final[.x, ],
      w    = 0,
      x    = opt$y_centers_final[.x, 1, ],
      y    = opt$y_centers_final[.x, 2, ],
      z    = opt$y_centers_final[.x, 3, ]
    )))
  })

  res <- list(
    qts_aligned = as_qts_sample(purrr::imap(x, ~ {
      .x$time <- opt$x_final[.y, ]
      .x
    })),
    qts_centers = centers,
    best_kma_result = opt
  )

  class(res) <- "kma_qts"
  res
}

#' QTS K-Means Visualization
#'
#' @param x An object of class `kma_qts` as produced by the [kmeans()] function.
#' @param ... Further arguments to be passed to other methods.
#'
#' @return The [plot.kma_qts()] method does not return anything while the
#'   [autoplot.kma_qts()] method returns a [ggplot2::ggplot] object.
#'
#' @importFrom graphics plot
#' @export
#'
#' @examples
#' res_kma <- kmeans(vespa64$igp, k = 2, nstart = 1)
#' plot(res_kma)
#' ggplot2::autoplot(res_kma)
plot.kma_qts <- function(x, ...) {
  print(autoplot(x, ...))
}

#' @importFrom ggplot2 autoplot .data
#' @export
#' @rdname plot.kma_qts
autoplot.kma_qts <- function(x, ...) {
  data <- as_qts_sample(c(x$qts_centers, x$qts_aligned))
  n <- length(x$qts_aligned)
  k <- length(x$qts_centers)
  memb <- c(1:k, x$best_kma_result$labels)
  high <- c(rep(TRUE, k), rep(FALSE, n))
  autoplot(data, memberships = memb, highlighted = high) +
    ggplot2::labs(
      title = "K-Means Alignment Clustering Results",
      subtitle = cli::pluralize("Using {k} cluster{?s}")
    )
}
