library(squat)

filter_by_pca <- function(qts_sample, component_ids, remove_components = FALSE) {
  N <- length(qts_sample)
  pca_res <- prcomp(vespa64$igp, M = max(component_ids))
  log_qts_sample <- purrr::map(1:N, ~ {
    tibble::tibble(
      time = qts_sample[[1]]$time,
      w = rep(0, 101),
      x = as.numeric(pca_res$tpca$scores[.x, component_ids] %*%
                       pca_res$tpca$functions[[1]]@X[component_ids, ] +
                       pca_res$tpca$meanFunction[[1]]@X),
      y = as.numeric(pca_res$tpca$scores[.x, component_ids] %*%
                       pca_res$tpca$functions[[2]]@X[component_ids, ] +
                       pca_res$tpca$meanFunction[[2]]@X),
      z = as.numeric(pca_res$tpca$scores[.x, component_ids] %*%
                       pca_res$tpca$functions[[3]]@X[component_ids, ] +
                       pca_res$tpca$meanFunction[[3]]@X)
    )
  })

  if (remove_components) {
    log_qts_sample <- purrr::map2(log(qts_sample), log_qts_sample, ~ {
      res <- .x
      res[, 2:5] <- res[, 2:5] - .y[, 2:5]
      res
    })
  }

  exp(as_qts_sample(log_qts_sample))
}

correctedQTS <- filter_by_pca(vespa64$igp, 1:2, remove_components = TRUE)

plot(vespa64$igp, memberships = vespa64$V)
plot(correctedQTS, memberships = vespa64$V)

plot(vespa64$igp, memberships = vespa64$P)
plot(correctedQTS, memberships = vespa64$P)

pca_res_corr <- prcomp(correctedQTS)
p <- autoplot(pca_res_corr, what = "scores")
p + geom_point(aes(color = vespa64$V))
p + geom_point(aes(color = vespa64$P))
