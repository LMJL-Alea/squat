#' Tangent PCA for QTS Data
#'
#' @param qts_list A list specifying the sample of QTS.
#' @param M An integer value specifying the number of principal component to
#'   compute. Defaults to `5L`.
#' @param fit A boolean specifying whether the resulting `qtsTPCA` object should
#'   store a reconstruction of the sample from the retained PCs. Defaults to
#'   `FALSE`.
#'
#' @return An object of class `qtsTPCA` which is a list with the following
#'   components:
#' - An object of class `MFPCAfit` as produced by the function
#' \code{\link[MFPCA]{MFPCA}},
#' - A \code{\link[tibble]{tibble}} containing the mean QTS,
#' - A list of \code{\link[tibble]{tibble}}s containing the required principal
#' components.
#'
#' @export
#'
#' @examples
#' pca_res <- tpca_qts(vespa$igp[1:8])
tpca_qts <- function(qts_list, M = 5, fit = FALSE) {
  check_common_grid <- qts_list |>
    purrr::map("time") |>
    purrr::reduce(rbind) |>
    apply(MARGIN = 2, FUN = stats::var) |>
    sum()
  if (check_common_grid > 0)
    cli::cli_abort("All input QTS should be evaluated on the same grid.")
  grid <- qts_list[[1]]$time
  qts_log <- purrr::map(qts_list, log_qts)
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
  uniExpansions <- purrr::map(1:3, ~ list(type = "uFPCA"))
  tpca <- MFPCA::MFPCA(mfd, M = M, uniExpansions = uniExpansions, fit = fit)
  mean_qts <- tpca$meanFunction |>
    purrr::map(~ as.numeric(.x@X)) |>
    purrr::set_names(c("x", "y", "z")) |>
    tibble::as_tibble()
  mean_qts$time <- grid
  mean_qts$w <- 0
  mean_qts <- exp_qts(mean_qts[c(4, 5, 1:3)])
  res <- list(
    tpca = tpca,
    mean_qts = mean_qts,
    principal_qts = tpca$functions |>
      purrr::map(~ purrr::array_tree(.x@X, margin = 1)) |>
      purrr::transpose() |>
      purrr::map(~ tibble::tibble(
        time = grid,
        w = 0,
        x = .x[[1]],
        y = .x[[2]],
        z = .x[[3]]
      )) |>
      purrr::map(exp_qts)
  )
  class(res) <- "qtsTPCA"
  res
}

#' @importFrom ggplot2 autoplot .data
autoplot.qtsTPCA <- function(x, what = "PC1", ...) {
  dots <- list(...)
  if (substr(what, 1, 2) == "PC") {
    component <- as.numeric(substr(what, 3, nchar(what)))
    if ("original_space" %in% names(dots))
      plot_tpca_component(x, component = component, original_space = dots$original_space)
    else {
      cli::cli_inform("The {.code original_space} boolean argument is not specified. Defaulting to {.field TRUE}.")
      plot_tpca_component(x, component = component, original_space = TRUE)
    }
  } else if (what == "scores") {
    if ("plane" %in% names(dots))
      plot_tpca_scores(x, plane = dots$plane)
    else {
      cli::cli_inform("The {.code plane} length-2 integer vector argument is not specified. Defaulting to {.field 1:2}.")
      plot_tpca_scores(x, plane = 1:2)
    }
  } else
    cli::cli_abort("The {.code what} argument should be either {.field scores} or a principal component specified starting with {.field PC}.")
}

#' Visualization of Tangent PCA for QTS
#'
#' @param x An object of class `qtsTPCA` as produced by \code{\link{tpca_qts}}.
#' @param what A string specifying what kind of visualization the user wants to
#'   perform. Choices are words starting with `PC` and ending with a PC number
#'   (in which case the mean QTS is displayed along with its perturbations due
#'   to the required PC) or `scores` (in which case individuals are projected on
#'   the required plance). Defaults to `PC1`.
#' @param ... If `what = "PC?"`, the user can specify whether to plot the QTS in
#'   the tangent space or in the original space by providing a boolean argument
#'   `original_space` which defaults to `TRUE`. If `what = "scores"`, the user
#'   can specify the plane onto which the individuals will be projected by
#'   providing a length-2 integer vector argument `plane` which defaults to
#'   `1:2`.
#'
#' @return NULL
#' @export
#'
#' @importFrom graphics plot
#'
#' @examples
#' pca_res <- tpca_qts(vespa$igp[1:8])
#' plot(pca_res)
plot.qtsTPCA <- function(x, what = "PC1", ...) {
  print(autoplot(x, what = what, ...))
}

plot_tpca_component <- function(tpca, component = 1, original_space = TRUE) {
  plot_mean <- log_qts(tpca$mean_qts)
  plot_cp <- log_qts(tpca$principal_qts[[component]])
  K <- stats::median(abs(tpca$tpca$scores[, component]))
  plot_lb <- plot_mean
  plot_lb[3:5] <- plot_lb[3:5] - K * plot_cp[3:5]
  plot_ub <- plot_mean
  plot_ub[3:5] <- plot_ub[3:5] + K * plot_cp[3:5]
  if (original_space) {
    plot_mean <- exp_qts(plot_mean)
    plot_lb <- exp_qts(plot_lb)
    plot_ub <- exp_qts(plot_ub)
  } else {
    plot_mean$w <- NULL
    plot_lb$w <- NULL
    plot_ub$w <- NULL
  }
  plot_mean$col <- "mean"
  plot_lb$col <- cli::pluralize("mean - med(|scores|) * PC{component}")
  plot_ub$col <- cli::pluralize("mean + med(|scores|) * PC{component}")
  rbind(plot_mean, plot_lb, plot_ub) |>
    tidyr::pivot_longer(-c(.data$time, .data$col)) |>
    ggplot2::ggplot(ggplot2::aes(x = .data$time, y = .data$value, color = .data$col)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(ggplot2::vars(.data$name), ncol = 1, scales = "free") +
    ggplot2::theme_linedraw() +
    ggplot2::labs(
      title = cli::pluralize("Mean QTS perturbed by PC{component}"),
      subtitle = cli::pluralize("Percentage of variance explained: {round(tpca$tpca$values[component] * 100, digits = 1)}%"),
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
  tibble::tibble(x = scores[, 1], y = scores[, 2]) |>
    ggplot2::ggplot(ggplot2::aes(.data$x, .data$y, label = 1:n)) +
    ggplot2::geom_point() +
    ggrepel::geom_label_repel() +
    ggplot2::theme_linedraw() +
    ggplot2::labs(
      title = cli::pluralize("Individuals projected on the PC{plane[1]}-{plane[2]} plane"),
      subtitle = cli::pluralize("Combined percentage of variance explained: {round(sum(tpca$tpca$values[plane]) * 100, digits = 1)}%"),
      x = cli::pluralize("PC{plane[1]} ({round(sum(tpca$tpca$values[plane[1]]) * 100, digits = 1)}%)"),
      y = cli::pluralize("PC{plane[2]} ({round(sum(tpca$tpca$values[plane[2]]) * 100, digits = 1)}%)")
    )
}
