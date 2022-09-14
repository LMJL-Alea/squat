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

qts_mean <- readRDS("data-raw/qts_mean.rds")
qts_mean <- purrr::map(qts_mean, squat::straighten_qts)
vespa <- purrr::map_df(names(qts_mean), parse_name)
vespa$igp <- `names<-`(qts_mean, NULL)

usethis::use_data(vespa, overwrite = TRUE, compress = "xz", version = 3)
