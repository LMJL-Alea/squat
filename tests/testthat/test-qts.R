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
  withr::with_seed(1234, {
    expect_snapshot(derivative_qts(vespa64$igp[[1]]))
  })
})

test_that("The function rnorm_qts() works", {
  withr::with_seed(1234, {
    expect_snapshot(rnorm_qts(1, vespa64$igp[[1]]))
  })
})

test_that("The function scale_qts() works (center = TRUE, by_row = FALSE, keep_summary_stats = FALSE)", {
  qts_list <- scale_qts(
    qts_list = vespa64$igp,
    center = TRUE,
    standardize = TRUE,
    by_row = FALSE,
    keep_summary_stats = FALSE
  )
  expect_snapshot(qts_list[[1]])
})

test_that("The function scale_qts() works (center = FALSE, by_row = FALSE, keep_summary_stats = FALSE)", {
  qts_list <- scale_qts(
    qts_list = vespa64$igp,
    center = FALSE,
    standardize = TRUE,
    by_row = FALSE,
    keep_summary_stats = FALSE
  )
  expect_equal(qts_list, vespa64$igp)
})

test_that("The function scale_qts() works (center = FALSE, by_row = FALSE, keep_summary_stats = TRUE)", {
  qts_list <- scale_qts(
    qts_list = vespa64$igp,
    center = FALSE,
    standardize = TRUE,
    by_row = FALSE,
    keep_summary_stats = TRUE
  )
  expect_equal(qts_list, list(
    qts_list = vespa64$igp,
    mean_values = NA,
    sd_values = NA
  ))
})

test_that("The function scale_qts() works (center = TRUE, by_row = TRUE, keep_summary_stats = FALSE)", {
  qts_list <- scale_qts(
    qts_list = vespa64$igp,
    center = TRUE,
    standardize = TRUE,
    by_row = TRUE,
    keep_summary_stats = FALSE
  )
  expect_snapshot(qts_list[[1]])
})
