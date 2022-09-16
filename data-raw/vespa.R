## code to prepare `vespa` dataset goes here

parse_name <- function(x) {
  dvalue <- substr(x, 2, 4)
  vvalue <- substr(x, 7, 7)
  evalue <- substr(x, 9, 9)
  rvalue <- substr(x, 11, 11)
  svalue <- substr(x, 13, 13)
  pvalue <- substr(x, 15, 15)
  avalue <- substr(x, 17, 17)
  list(
    "V" = vvalue,
    "E" = evalue,
    "S" = svalue,
    "P" = pvalue,
    "A" = avalue,
    "R" = rvalue
  )
}

qts_mean <- readRDS("data-raw/qts_mean.rds") |>
  purrr::map(squat::as_qts) |>
  purrr::map(squat::straighten_qts)
vespa <- purrr::map_df(names(qts_mean), parse_name)
vespa$igp <- `names<-`(qts_mean, NULL)

usethis::use_data(vespa, overwrite = TRUE, compress = "xz", version = 3)

vespa64 <- vespa |>
  mutate(q = map(igp, ~ pmap(list(.x$w, .x$x, .x$y, .x$z), c))) |>
  select(-igp) |>
  unnest(q) |>
  group_by(V, E, S, P, A, R) |>
  mutate(id = 1:n()) |>
  group_by(V, E, S, P, A, id) |>
  summarise(q = list(squat:::gmean(q))) |>
  group_by(V, E, S, P, A) |>
  summarise(igp = list(tibble(
    time = 0:100,
    w = map_dbl(q, 1),
    x = map_dbl(q, 2),
    y = map_dbl(q, 3),
    z = map_dbl(q, 4)
  ))) |>
  ungroup() |>
  mutate(igp = map(igp, as_qts))

usethis::use_data(vespa64, overwrite = TRUE, compress = "xz", version = 3)
