test_that("the function qts2dts() works", {
  expect_snapshot(qts2dts(vespa64$igp[[1]], vespa64$igp[[2]]))
})
