#' Grid regularization
#'
#' This function makes sure that a quaternion time series is evaluated on a grid
#' with fixed step size.
#'
#' @param x A numeric vector providing the original evaluation grid.
#' @param y A \code{4 x length(x)} matrix providing the original QTS.
#' @param xmin A real scalar specifying the minimum of the output evaluation
#'   grid (default: minimum of input evaluation grid).
#' @param xmax A real scalar specifying the maximum of the output evaluation
#'   grid (default: maximum of input evaluation grid).
#' @param out_size An integer specifying the size of the output evaluation grid.
#'   By default, it takes the same size as the input evaluation grid.
#'
#' @return A \code{4 x length(x)} matrix providing the regularized QTS.
#' @export
#'
#' @examples
regularize_grid <- function(x, y, xmin = min(x), xmax = max(x), out_size = 0) {
  if (xmin < min(x) | xmax > max(x))
    stop("It is not possible to interpolate outside the evaluation grid range.")
  RegularizeGrid(x, y, xmin, xmax, out_size)
}
