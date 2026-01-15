test_that("the function qts2dts() works", {
  expect_snapshot(qts2dts(vespa64$igp[[1]], vespa64$igp[[2]]))
})

test_that("the function qts2nts() works", {
  expect_snapshot(qts2nts(vespa64$igp[[1]], disable_normalization = FALSE))
  expect_snapshot(qts2nts(vespa64$igp[[1]], disable_normalization = TRUE))
})

test_that("the function qts2ats() works", {
  expect_snapshot(qts2ats(vespa64$igp[[1]], disable_normalization = FALSE))
  expect_snapshot(qts2ats(vespa64$igp[[1]], disable_normalization = TRUE))
})

test_that("the function qts2avvts() works", {
  expect_snapshot(qts2avvts(vespa64$igp[[1]]))
  expect_snapshot(qts2avmts(vespa64$igp[[1]]))
})

test_that("the function qts2aavts() works", {
  expect_snapshot(qts2aavts(vespa64$igp[[1]]))
  expect_snapshot(qts2aamts(vespa64$igp[[1]]))
})

test_that("the function qts2sqts() works", {
  expect_snapshot(qts2sqts(vespa64$igp[[1]]))
})
