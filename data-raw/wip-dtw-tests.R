library(future)
system.time(
  D1 <- dist(vespa64$igp)
)
plan(multisession, workers = 8)
system.time(
  D2 <- dist(vespa64$igp)
)
plan(sequential)
waldo::compare(D1, D2)

DTW(vespa64$igp[[1]], vespa64$igp[[2]])$distance
DTWi(vespa64$igp[[1]], vespa64$igp[[2]])$distance
DTWi(vespa64$igp[[2]], vespa64$igp[[1]])$distance

# DTW on mean-centered
D1 <- dist(vespa64$igp, metric = "dtw")
plan(multisession, workers = 8)
# DTWi on mean-centered
D2 <- dist(vespa64$igp, metric = "dtw")
# DTWi on mean-centered + sbplx
D5 <- dist(vespa64$igp, metric = "dtw")
# DTWi on mean-centered + k=2
D6 <- dist(vespa64$igp, metric = "dtw")
plan(sequential)
barplot(as.numeric(D2 - D6))
# DTW on first-centered
D3 <- dist(reorient(vespa64$igp), metric = "dtw")
plan(multisession, workers = 8)
# DTWi on first-centered
D4 <- dist(reorient(vespa64$igp), metric = "dtw")
future::plan(future::sequential)
barplot(as.numeric(D3 - D4))
barplot(as.numeric(D2 - D4))
summary(as.numeric(D2 - D4))
boxplot(as.numeric(D2 - D4))
