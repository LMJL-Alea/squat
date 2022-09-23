test_that("Function distDTW() works", {
  D <- distDTW(vespa64$igp)
  expect_snapshot(D)
})
