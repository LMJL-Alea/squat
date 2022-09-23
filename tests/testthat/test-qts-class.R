test_that("Functions related to the QTS class work", {
  qts1 <- vespa64$igp[[1]]
  expect_true(is_qts(qts1))
  qts2 <- as_qts(qts1)
  expect_true(is_qts(qts2))
  expect_equal(qts1, qts2)
  qts3 <- qts1
  class(qts3) <- class(qts3)[-1]
  expect_false(is_qts(qts3))
  qts3 <- as_qts(qts1)
  expect_equal(qts1, qts3)
})

test_that("The function derivative_qts() works", {
  expect_snapshot(derivative_qts(vespa64$igp[[1]]))
})

test_that("Logarithm and exponential for QTS work", {
  x <- log_qts(vespa64$igp[[1]])
  expect_snapshot(x)
  y <- exp_qts(x)
  expect_equal(y, vespa64$igp[[1]])
})

test_that("Function reorient_qts() works", {
  expect_snapshot(reorient_qts(vespa64$igp[[1]], disable_normalization = FALSE))
  expect_snapshot(reorient_qts(vespa64$igp[[1]], disable_normalization = TRUE))
})

test_that("Function normalize_qts() works", {
  expect_snapshot(normalize_qts(vespa64$igp[[1]]))
})

test_that("Function centring_qts() works (standardize = FALSE, keep_summary_stats = FALSE)", {
  expect_snapshot(centring_qts(
    x = vespa64$igp[[1]],
    standardize = FALSE,
    keep_summary_stats = FALSE
  ))
})

test_that("Function centring_qts() works (standardize = TRUE, keep_summary_stats = FALSE)", {
  expect_snapshot(centring_qts(
    x = vespa64$igp[[1]],
    standardize = TRUE,
    keep_summary_stats = FALSE
  ))
})

test_that("Function centring_qts() works (standardize = FALSE, keep_summary_stats = TRUE)", {
  expect_snapshot(centring_qts(
    x = vespa64$igp[[1]],
    standardize = FALSE,
    keep_summary_stats = TRUE
  ))
})

test_that("Function resample_qts() works", {
  expect_snapshot(resample_qts(vespa64$igp[[1]]))
})

test_that("Function smooth_qts() works", {
  expect_snapshot(smooth_qts(vespa64$igp[[1]]))
})

test_that("Visualization code for QTS work", {
  p <- ggplot2::autoplot(vespa64$igp[[1]])
  expect_equal(dim(p$data), c(404, 3))
  p <- ggplot2::autoplot(vespa64$igp[[1]], change_points = c(10, 80))
  expect_equal(dim(p$data), c(404, 3))
})

test_that("Visualization functions for QTS work", {
  skip_if_not_installed("vdiffr")
  skip_on_covr()
  skip_on_ci()
  vdiffr::expect_doppelganger(
    title = "QTS plot",
    fig = plot(vespa64$igp[[1]])
  )
  vdiffr::expect_doppelganger(
    title = "QTS plot with change points",
    fig = plot(vespa64$igp[[1]], change_points = c(10, 80))
  )
})
