# socialroulette
R package for partitioning individuals into groups
```
# Distribute 5 ppl into groups of at least 2
frame <- tibble(id=sprintf("id%.3d@math.su.se", 1:5), date=Sys.Date())
rsocialroulette(frame, past_sessions=NULL, m=2, algorithm="srs")


# Create a history of previous sessions as well as the current frame from a population of 100
ppl <- tibble(id=sprintf("id%.3d@math.su.se", 1:100))
session_dates <- seq(as.Date("2021-04-01"), length.out=4, by="1 week")

# Simulate changing participation each week (with prob 0.7 they will attend)
frames <- map_df( session_dates, ~ ppl %>% slice_sample(n = rbinom(1,nrow(.), prob=0.7)) %>% mutate(date=.x))

# Frame for the current session is the last one
current_frame <- frames %>% filter(date == max(date))

# Simulate the partitions for the last 3 weeks
past_sessions <- frames %>% filter(date < max(session_dates)) %>%
  group_split(date) %>%
  map(~rsocialroulette(current_frame=.x, past_sessions=NULL, m=4)) %>%
  setNames(session_dates[1:3])
```

## Make the partition, which tries to bring together those who have not met in a long time
```
rsocialroulette(current_frame, past_session, m=4, algorithm="mdgp")
```
