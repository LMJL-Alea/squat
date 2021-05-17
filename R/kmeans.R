#' QTS K-Mean Alignment Algorithm
#'
#' This function massages the input quaternion time series to feed them into the
#' k-mean alignment algorithm for jointly clustering and aligning the input QTS.
#'
#' @param qts_list A list of QTS stored as \code{\link[tibble]{tibble}}s with
#'   columns `time`, `w`, `x`, `y` and `z`.
#' @param k An integer specifying the number of clusters to be formed. Defaults
#'   to `1L`.
#' @param centroid A string specifying which type of centroid should be used.
#'   Choices are `mean` and `medoid`. Defaults to `mean`.
#' @param diss A string specifying which type of dissimilarity should be used.
#'   Choices are `l1` and `pearson`. Defaults to `l2`.
#' @param warping A string specifying which class of warping functions should be
#'   used. Choices are `none`, `shift`, `dilation` and `affine`. Defaults to
#'   `affine`.
#' @param iter_max An integer specifying the maximum number of iterations for
#'   terminating the k-mean algorithm. Defautls to `20L`.
#' @param nstart An integer specifying the number of random restart for making
#'   the k-mean results more robust. Defaults to `1000L`.
#' @param ncores An integer specifying the number of cores to run the multiple
#'   restarts of the k-mean algorithm in parallel. Defaults to `1L`.
#'
#' @return A \code{\link[fdakmapp]{kma}} object storing the results of the best
#'   k-mean alignment algorithm run.
#' @export
#'
#' @examples
#' # TO DO
kmeans_qts <-function(qts_list,
                      k = 1,
                      centroid = "mean",
                      diss = "l2",
                      warping = "affine",
                      iter_max = 20,
                      nstart = 1000,
                      ncores = 1L) {
  q <- lapply(qts_list, squat::log_qts)
  q <- lapply(q, function(.qts) {
      rbind(.qts$x, .qts$y, .qts$z)
    })
  t <- lapply(qts_list, function(.qts) .qts$time)

  # Prep data
  n <- length(q)
  d <- dim(q[[1]])[1]
  p <- dim(q[[1]])[2]

  if (is.null(t))
    x <- 0:(p-1)
  else
    x <- matrix(nrow = n, ncol = p)

  y <- array(dim = c(n, d, p))
  for (i in 1:n) {
    y[i, , ] <- q[[i]]
    if (!is.null(t)) {
      x[i, ] <- t[[i]]
    }
  }

  cl <- NULL
  if (ncores > 1L)
    cl <- parallel::makeCluster(ncores)

  B <- choose(n, k)

  if (nstart > B)
    init <- utils::combn(n, k, simplify = FALSE)
  else
    init <- replicate(nstart, sample.int(n, k), simplify = FALSE)

  solutions <- pbapply::pblapply(
    X = init,
    FUN = function(.init) {
      fdakmapp::kma(
        x,
        y,
        seeds = .init,
        n_clust = k,
        center_method = centroid,
        warping_method = warping,
        dissimilarity_method = diss,
        use_verbose = FALSE,
        space = 0,
        warping_options = c(0.1, 0.1),
        use_fence = FALSE
      )
    },
    cl = cl
  )

  if (!is.null(cl))
    parallel::stopCluster(cl)

  wss_vector <- sapply(solutions, function(.sol) {
    sum(.sol$final_dissimilarity)
  })

  opt <- solutions[[which.min(wss_vector)]]

  qts_center <- tibble(
    time = opt$x_centers_final[1, ],
    w    = 0,
    x    = opt$y_centers_final[1, 1, ],
    y    = opt$y_centers_final[1, 2, ],
    z    = opt$y_centers_final[1, 3, ]
  )

  list(
    qts_aligned = lapply(1:length(qts_list), function(.idx) {
      qts_list[[.idx]]$time <- opt$x_final[.idx, ]
      qts_list[[.idx]]
    }),
    qts_center = exp_qts(qts_center),
    best_kma_result = opt
  )
}
