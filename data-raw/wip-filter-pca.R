library(squat)

plot(vespa64$igp, memberships = vespa64$V)
plot(vespa64$igp, memberships = vespa64$P)

pca_res <- prcomp(vespa64$igp)
p <- autoplot(pca_res, what = "scores")
p + geom_point(color = vespa64$V)
p + geom_point(color = vespa64$P)

filter_by_pca <- function(qts_sample, component_ids, remove_components = FALSE) {
  N <- length(qts_sample)
  pca_res <- prcomp(qts_sample, M = max(component_ids))
  log_qts_sample <- purrr::map(1:N, ~ {
    df <- tibble::tibble(
      time = qts_sample[[1]]$time,
      w = rep(0, 101),
      x = as.numeric(pca_res$tpca$scores[.x, component_ids] %*%
                       pca_res$tpca$functions[[1]]@X[component_ids, ]),
      y = as.numeric(pca_res$tpca$scores[.x, component_ids] %*%
                       pca_res$tpca$functions[[2]]@X[component_ids, ]),
      z = as.numeric(pca_res$tpca$scores[.x, component_ids] %*%
                       pca_res$tpca$functions[[3]]@X[component_ids, ])
    )
    if (!remove_components) {
      df$x <- df$x + pca_res$tpca$meanFunction[[1]]@X
      df$y <- df$y + pca_res$tpca$meanFunction[[2]]@X
      df$z <- df$z + pca_res$tpca$meanFunction[[3]]@X
    }
    df
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

correctedQTS <- filter_by_pca(vespa64$igp, 2, remove_components = TRUE)

D1 <- distKMA(vespa64$igp)
D2 <- distKMA(correctedQTS)

mod1 <- vegan::adonis2(D1 ~ V+E+S+P, data = vespa64, add = "lingoes")
mod2 <- vegan::adonis2(D2 ~ V+E+S+P, data = vespa64, add = "lingoes")

mod1
mod2

plot(vespa64$igp, memberships = vespa64$V)
plot(vespa64$igp, memberships = vespa64$P)
plot(correctedQTS, memberships = vespa64$V)

plot(vespa64$igp, memberships = vespa64$P)
plot(correctedQTS, memberships = vespa64$P)

pca_res_corr <- prcomp(correctedQTS)
p <- autoplot(pca_res_corr, what = "scores")
p + geom_point(aes(color = vespa64$V))
p + geom_point(aes(color = vespa64$P))
p + geom_point(aes(color = vespa64$S))

vespa_corr <- vespa64
vespa_corr$igp <- correctedQTS
vespa_corr <- vespa_corr |>
  group_by(V, E, S, A) |>
  summarise(igp = list(mean(igp))) |>
  ungroup() |>
  mutate(igp = as_qts_sample(igp))

pca_res_corr <- prcomp(vespa_corr$igp)
p <- autoplot(pca_res_corr, what = "scores")
p + geom_point(aes(color = vespa_corr$V))

correctedQTS <- filter_by_pca(vespa_corr$igp, 1, remove_components = TRUE)

D3 <- dist(correctedQTS, metric = "dtw", warping_class = "none")

mod3 <- vegan::adonis2(D3 ~ V+E+S, data = vespa_corr, add = "lingoes")
mod3
vegan::adonis2(D3 ~ E, data = vespa_corr, add = "lingoes")


plot(vespa64$igp, memberships = vespa64$V)
plot(correctedQTS, memberships = vespa_corr$V)

p_position <- autoplot(pca_res, what = "PC2")
p_position + labs(title = "Sensor Position Modulation")
p_volunteer <- autoplot(pca_res_corr, what = "PC1")
p_volunteer + labs(title = "Natural Gait Modulation")

Xorig <- map(1:3, ~ {
  vespa64$igp |>
    exp() |>
    map(function(qts) qts[[2+.x]]) |>
    do.call(what = cbind, args = _)
})
fdANOVA::fmanova.ptbfr(Xorig, group.label = vespa64$V, int = c(0, 100))

Xpartial <- map(1:3, ~ {
  vespa_corr$igp |>
    exp() |>
    map(function(qts) qts[[2+.x]]) |>
    do.call(what = cbind, args = _)
})
fdANOVA::fmanova.ptbfr(Xpartial, group.label = vespa_corr$V, int = c(0, 100))

Xcorr <- map(1:3, ~ {
  correctedQTS |>
    exp() |>
    map(function(qts) qts[[2+.x]]) |>
    do.call(what = cbind, args = _)
})
fdANOVA::fmanova.ptbfr(Xcorr, group.label = vespa_corr$V, int = c(0, 100))
