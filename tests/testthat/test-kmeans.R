test_that("multiplication works", {
  res_kma <- kmeans_qts(vespa64$igp, k = 2, nstart = 1)
  expect_equal(res_kma$best_kma_result$n_clust, 2)
})
