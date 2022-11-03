test_that("multiplication works", {
  withr::with_seed(1234, {
    res_kma <- kmeans(vespa64$igp, k = 2, nstart = 1)
  })
  expect_equal(res_kma$best_kma_result$n_clust, 2)
})

test_that("Visualization code for k-means work", {
  withr::with_seed(1234, {
    res_kma <- kmeans(vespa64$igp, k = 2, nstart = 1)
  })
  p <- ggplot2::autoplot(res_kma)
  expect_equal(dim(p$data), c(808, 6))
})

test_that("Visualization functions for PCA work", {
  skip_if_not_installed("vdiffr")
  skip_on_covr()
  skip_on_ci()
  withr::with_seed(1234, {
    res_kma <- kmeans(vespa64$igp, k = 2, nstart = 1)
  })
  vdiffr::expect_doppelganger(
    title = "K-means plot",
    fig = plot(res_kma)
  )
})
