test_that("The function tpca_qts() works", {
  res_pca <- tpca_qts(vespa64$igp)
  expect_snapshot(res_pca)
})

test_that("Visualization functions for PCA work", {
  skip_if_not_installed("vdiffr")
  skip_on_covr()
  skip_on_ci()
  res_pca <- tpca_qts(vespa64$igp)
  vdiffr::expect_doppelganger(
    title = "PC plot",
    fig = plot(res_pca, what = "PC1")
  )
  vdiffr::expect_doppelganger(
    title = "Score plot",
    fig = plot(res_pca, what = "scores")
  )
  p <- ggplot2::autoplot(res_pca, what = "scores")
  vdiffr::expect_doppelganger(
    title = "Colored score plot",
    fig = p + ggplot2::geom_point(ggplot2::aes(color = vespa64$V))
  )
})
