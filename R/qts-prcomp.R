#' PCA for QTS Sample
#'
#' @param x An object of class [qts_sample].
#' @param M An integer value specifying the number of principal component to
#'   compute. Defaults to `5L`.
#' @param fit A boolean specifying whether the resulting `prcomp_qts` object
#'   should store a reconstruction of the sample from the retained PCs. Defaults
#'   to `FALSE`.
#' @param ... Arguments passed to or from other methods.
#'
#' @return An object of class `prcomp_qts` which is a list with the following
#'   components:
#' - `tpca`: An object of class `MFPCAfit` as produced by the function
#' [MFPCA::MFPCA()],
#' - `var_props`: A numeric vector storing the percentage of variance explained
#' by each PC,
#' - `mean_qts`: An object of class [qts] containing the mean QTS,
#' - `principal_qts`: A list of [qts]s containing the required principal
#' components.
#'
#' @importFrom stats prcomp
#' @export
#'
#' @examples
#' res_pca <- prcomp(vespa64$igp[1:16])
prcomp.qts_sample <- function(x, M = 5, fit = FALSE, ...) {
  check_common_grid <- x |>
    purrr::map("time") |>
    purrr::reduce(rbind) |>
    apply(MARGIN = 2, FUN = stats::var) |>
    sum()
  if (check_common_grid > 0)
    cli::cli_abort("All input QTS should be evaluated on the same grid.")
  grid <- x[[1]]$time
  npc_max <- min(length(x), length(grid)) - 1
  if (M > npc_max)
    cli::cli_abort("The maximum number of principal component is {npc_max}. Please choose a value of {.arg M} smaller or equal to that value.")
  qts_log <- log(x)
  fd_x <- funData::funData(
    argvals = grid,
    X = qts_log |>
      purrr::map("x") |>
      purrr::reduce(rbind)
  )
  fd_y <- funData::funData(
    argvals = grid,
    X = qts_log |>
      purrr::map("y") |>
      purrr::reduce(rbind)
  )
  fd_z <- funData::funData(
    argvals = grid,
    X = qts_log |>
      purrr::map("z") |>
      purrr::reduce(rbind)
  )
  mfd <- funData::multiFunData(fd_x, fd_y, fd_z)
  uniExpansions <- purrr::map(1:3, ~ list(type = "splines1Dpen", k = npc_max))
  tpca <- MFPCA::MFPCA(mfd, M = M, uniExpansions = uniExpansions, fit = fit)
  # uniExpansions <- purrr::map(1:3, ~ list(type = "splines1D", k = M))
  # tpca <- MFPCA::MFPCA(mfd, M = M, uniExpansions = uniExpansions, fit = fit, approx.eigen = TRUE)
  if (M == npc_max)
    tot_var <- sum(tpca$values)
  else {
    # uniExpansions <- purrr::map(1:3, ~ list(type = "splines1Dpen", k = npc_max))
    tpca_full <- MFPCA::MFPCA(mfd, M = npc_max, uniExpansions = uniExpansions, fit = FALSE)
    # uniExpansions <- purrr::map(1:3, ~ list(type = "splines1D", k = npc_max))
    # tpca_full <- MFPCA::MFPCA(mfd, M = npc_max, uniExpansions = uniExpansions, fit = FALSE, approx.eigen = TRUE)
    tot_var <- sum(tpca_full$values)
  }
  mean_qts <- tpca$meanFunction |>
    purrr::map(~ as.numeric(.x@X)) |>
    purrr::set_names(c("x", "y", "z")) |>
    tibble::as_tibble()
  mean_qts$time <- grid
  mean_qts$w <- 0
  mean_qts <- exp(as_qts(mean_qts[c(4, 5, 1:3)]))
  res <- list(
    tpca = tpca,
    var_props = tpca$values / tot_var,
    mean_qts = mean_qts,
    principal_qts = tpca$functions |>
      purrr::map(~ purrr::array_tree(.x@X, margin = 1)) |>
      purrr::transpose() |>
      purrr::map(~ as_qts(tibble::tibble(
        time = grid,
        w = 0,
        x = .x[[1]],
        y = .x[[2]],
        z = .x[[3]]
      ))) |>
      purrr::map(exp)
  )
  class(res) <- "prcomp_qts"
  res
}

#' Plot for `prcomp_qts` objects
#'
#' This function creates a visualization of the results of the PCA applied on a
#' sample of QTS and returns the corresponding [ggplot2::ggplot] object which
#' enable further customization of the plot.
#'
#' @param object An object of class `prcomp_qts` as produced by the
#'   [prcomp.qts_sample()] method.
#' @param what A string specifying what kind of visualization the user wants to
#'   perform. Choices are words starting with `PC` and ending with a PC number
#'   (in which case the mean QTS is displayed along with its perturbations due
#'   to the required PC) or `scores` (in which case individuals are projected on
#'   the required plane). Defaults to `PC1`.
#' @param ... If `what = "PC?"`, the user can specify whether to plot the QTS in
#'   the tangent space or in the original space by providing a boolean argument
#'   `original_space` which defaults to `TRUE`. If `what = "scores"`, the user
#'   can specify the plane onto which the individuals will be projected by
#'   providing a length-2 integer vector argument `plane` which defaults to
#'   `1:2`.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @importFrom ggplot2 autoplot .data
#' @export
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' df <- as_qts_sample(vespa64$igp[1:16])
#' res_pca <- prcomp(df)
#'
#' # Plot the data points in a PC plane
#' # And color points according to a categorical variable
#' p <- ggplot2::autoplot(res_pca, what = "scores")
#' p + ggplot2::geom_point(ggplot2::aes(color = vespa64$V[1:16]))
autoplot.prcomp_qts <- function(object, what = "PC1", ...) {
  dots <- list(...)
  if (substr(what, 1, 2) == "PC") {
    component <- as.numeric(substr(what, 3, nchar(what)))
    if ("original_space" %in% names(dots))
      plot_tpca_component(object, component = component, original_space = dots$original_space)
    else {
      cli::cli_inform("The {.code original_space} boolean argument is not specified. Defaulting to {.field TRUE}.")
      plot_tpca_component(object, component = component, original_space = TRUE)
    }
  } else if (what == "scores") {
    if ("plane" %in% names(dots))
      plot_tpca_scores(object, plane = dots$plane)
    else {
      cli::cli_inform("The {.code plane} length-2 integer vector argument is not specified. Defaulting to {.field 1:2}.")
      plot_tpca_scores(object, plane = 1:2)
    }
  } else if (what == "variance") {
    screeplot(object)
  } else
    cli::cli_abort("The {.arg what} argument should be either {.field scores} or {.field variance} or a principal component specified starting with {.field PC}.")
}

#' Plot for `prcomp_qts` objects
#'
#' This function creates a visualization of the results of the PCA applied on a
#' sample of QTS **without** returning the plot data as an object.
#'
#' @param x An object of class `prcomp_qts` as produced by the
#'   [prcomp.qts_sample()] method.
#' @inheritParams autoplot.prcomp_qts
#'
#' @return NULL
#'
#' @importFrom graphics plot
#' @export
#' @examples
#' df <- as_qts_sample(vespa64$igp[1:16])
#' res_pca <- prcomp(df)
#'
#' # You can plot the effect of a PC on the mean
#' plot(res_pca, what = "PC1")
#'
#' # You can plot the data points in a PC plane
#' plot(res_pca, what = "scores")
plot.prcomp_qts <- function(x, what = "PC1", ...) {
  print(autoplot(x, what = what, ...))
}

#' @importFrom stats screeplot
#' @export
#' @rdname plot.prcomp_qts
screeplot.prcomp_qts <- function(x, ...) {
  plot_data <- tibble::tibble(
    lambda = x$var_props,
    m = seq_along(.data$lambda)
  )

  plot_data |>
    ggplot2::ggplot(ggplot2::aes(.data$m, .data$lambda)) +
    ggplot2::geom_col() +
    ggplot2::labs(
      title = "Screeplot of QTS PCA",
      x = "Principal Component Index",
      y = "Percentage of variance explained"
    ) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::theme_linedraw()
}

plot_tpca_component <- function(tpca, component = 1, original_space = TRUE) {
  plot_mean <- log(tpca$mean_qts)
  plot_cp <- log(tpca$principal_qts[[component]])
  K <- stats::median(abs(tpca$tpca$scores[, component]))
  plot_lb <- plot_mean
  plot_lb[3:5] <- plot_lb[3:5] - K * plot_cp[3:5]
  plot_ub <- plot_mean
  plot_ub[3:5] <- plot_ub[3:5] + K * plot_cp[3:5]
  if (original_space) {
    plot_mean <- exp(plot_mean)
    plot_lb <- exp(plot_lb)
    plot_ub <- exp(plot_ub)
  } else {
    plot_mean$w <- NULL
    plot_lb$w <- NULL
    plot_ub$w <- NULL
  }
  plot_mean$col <- "mean"
  plot_lb$col <- cli::pluralize("mean - med(|scores|) * PC{component}")
  plot_ub$col <- cli::pluralize("mean + med(|scores|) * PC{component}")

  plot_data <- rbind(plot_mean, plot_lb, plot_ub) |>
    tidyr::pivot_longer(-c("time", "col"))

  plot_data |>
    ggplot2::ggplot(ggplot2::aes(x = .data$time, y = .data$value, color = .data$col)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(ggplot2::vars(.data$name), ncol = 1, scales = "free") +
    ggplot2::theme_linedraw() +
    ggplot2::labs(
      title = cli::pluralize("Mean QTS perturbed by PC{component}"),
      subtitle = cli::pluralize("Percentage of variance explained: {round(tpca$var_props[component] * 100, digits = 1)}%"),
      x = "Time (%)",
      y = "",
      color = ""
    )
}

plot_tpca_scores <- function(tpca, plane = 1:2) {
  if (length(plane) != 2)
    cli::cli_abort("The {.code plane} argument should be of length two.")
  scores <- tpca$tpca$scores[, plane]
  n <- nrow(scores)

  plot_data <- tibble::tibble(x = scores[, 1], y = scores[, 2])

  plot_data |>
    ggplot2::ggplot(ggplot2::aes(.data$x, .data$y, label = 1:n)) +
    ggplot2::geom_point() +
    ggrepel::geom_label_repel(seed = 1234) +
    ggplot2::theme_linedraw() +
    ggplot2::labs(
      title = cli::pluralize("Individuals projected on the PC{plane[1]}-{plane[2]} plane"),
      subtitle = cli::pluralize("Combined percentage of variance explained: {round(sum(tpca$var_props[plane]) * 100, digits = 1)}%"),
      x = cli::pluralize("PC{plane[1]} ({round(sum(tpca$var_props[plane[1]]) * 100, digits = 1)}%)"),
      y = cli::pluralize("PC{plane[2]} ({round(sum(tpca$var_props[plane[2]]) * 100, digits = 1)}%)")
    )
}
