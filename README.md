# socialroulette: A package to generate social groupings

A package for partitioning individuals into groups of a pre-specified size. This can be as simple as using simple random sampling (srs) to divide n individuals into groups of size at least m. If one keeps track of the past partitions then an additional
aim can be to try to maximize the time that people have a reunion, i.e. end up in the same group. This boils
down to an instance of the maximally diverse grouping problem.
See the [package website](https://hoehleatsu.github.io/socialroulette/)  for more information, documentation and examples. The package source code is available from [GitHub](https://github.com/hoehleatsu/socialroulette/).

```
# Install package
devtools::install_github("hoehleatsu/socialroulette")

# Distribute 5 ppl into groups of at least 2
frame <- tibble(id=sprintf("id%.3d@math.su.se", 1:5), date=Sys.Date())
rsocialroulette(frame, past_sessions=NULL, m=2, algorithm="srs")
```

For more information see the `get-started` vignette.
