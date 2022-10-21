sd_data <- vespa |>
  left_join(vespa64 |> rename(mean_igp = igp)) |>
  mutate(dist = map2(igp, mean_igp, DTW) |> map_dbl("normalizedDistance")) |>
  group_by(V, E, S, P, A) |>
  summarise(sd = sqrt(mean(dist^2))) |>
  filter(S %in% c("1", "3"))
sd_data |>
  ggplot(aes(S, sd, fill = V)) +
  geom_col(position = position_dodge()) +
  facet_grid(rows = vars(E), cols = vars(P))
