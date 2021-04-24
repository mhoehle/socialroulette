#' @import tidyverse
NULL

#' Maximally diverse grouping problem solver based on Lai and Hao (2016)
#'
#' The following
#' @param mdg_format_file Path to the file containing the MDGP specification in mdgplib format
#' @param time_limit Number of seconds to iteratively optimize each run
#' @return File name of the solution file
#'
#' @description The code calls C++ code by Xiangjing Lai and Jin-Kao Hao,
#' for this several temporary files are generated using `tempfile()`.
#' @details The maximally diverse grouping problem is about partitioning $n$ individuals
#' into groups of size at least $m$ while maximizing a sum of utility values computed
#' by summing the utility $u(i,j)$ over all individuals $i,j$ partitioned into the same group.
#' and summing this quantity over all groups.
#'
#' @author Xiangjing Lai and Jin-Kao Hao, R interface by M. Höhle
#' @useDynLib socialroulette mdgp
mdgp_solver <- function(mdgp_format_file, time_limit= 30) {
  #Sanity checks
  stopifnot(is.character(mdgp_format_file))
  stopifnot(file.exists(mdgp_format_file))

  #Temporary file for the output
  output_file <- tempfile()
  solution_file <- tempfile()

  #Call C++ code
  .C("mdgp",
     mdgp_format_file, #"lunch2.txt",
     output_file,
     solution_file,
     as.numeric(time_limit))

  #Read result into R compatible structure
  return(list(output_file=output_file, solution_file=solution_file))
}

#' Convert session specification to Lai and Hao (2016) MDGP solver input specification
#'
#' @param current_frame A pair-distance tibble. The frame needs to contain a column denoted `id`
#' @param past_sessions A list of sessions. Each sessions consists of a session specification in list format.
#' @param m The minimum group size, i.e. all groups have size at least m
#' @return File name of the generated specification file
#' @examples
#'  ppl <- tibble(id=sprintf("id%.3d", 1:100))
#'  session_dates <- seq(as.Date("2021-04-01"), length.out=4, by="1 week")
#'  frames <- map_df( session_dates, ~ ppl %>% slice_sample(n = rbinom(1,nrow(ppl), prob=0.7)) %>% mutate(date=.x))
#'  past_sessions <- frames %>% group_split(date) %>% map(~sample_groups(.x, m=4)) %>% setNames(session_dates) %>% .[1:3]
#'  current_frame <- frames %>% filter(date == max(date))
#'  spec_file <- write_mdgp_specfile(current_frame, past_sessions, m=4)
#'  mdgp_solver(spec_file)

write_mdgp_specfile <- function(current_frame, past_sessions, m) {
  # Check that the current_frame has an id column
  stopifnot( "id" %in% colnames(current_frame))

  # No. individuals in the frame
  n <- nrow(current_frame)
  # Number of resulting groups
  n_g <-  n %/% m
  # Group to assign each entry in frame into
  g <- rep(seq_len(n %/% m), length.out=n)
  # How many in each group? Note: can be way beyond m, even if m was selected
  n_groupsize <- table(table(g))
  # Vector containing min and max of each group
  g_spec <- sort(rep(as.numeric(names(n_groupsize)), times=n_groupsize), decreasing=TRUE)
  g_spec <- rep(g_spec, each=2)
  spec <- str_c(n, " ",n_g, " ds ", str_c(g_spec, collapse=" "))

  # Compute distances to last meet for each id in current_frame
  dist_frame <- sessions_to_distance(current_frame, past_sessions)

  # Make alternative index number (starting from zero) for MDGP-solver instead of the id
  ids <- current_frame %>% mutate(idx = row_number() - 1) %>% select(id, idx)
  # and add this to the dist_frame
  dist_frame <-  left_join(dist_frame, ids, by=c("id1"="id")) %>%
    left_join(ids, by=c("id2"="id")) %>%
    rename(idx.id1 = idx.x, idx.id2=idx.y) %>%
    #Make sure idx.id1 < idx.id2 as this is needed for the solver
    mutate( idx.id1.ord = pmin(idx.id1, idx.id2),
            idx.id2.ord = pmax(idx.id1, idx.id2)) %>%
    arrange(idx.id1.ord, idx.id2.ord)

  #Sanity check
  stopifnot( all(dist_frame$idx.id1.ord < dist_frame$idx.id2.ord))

  #Make a temporary file
  tmp_file <- str_c(tempfile(), ".txt")
  writeLines(spec, tmp_file)
  write_delim(dist_frame %>% select(idx.id1.ord,idx.id2.ord,dist), path=tmp_file, col_names=FALSE, append=TRUE)

  #Done - return filename
  return(tmp_file)
}

#' Read output from the Lai and Hao (2016) MDGP solver
#' @param file_name Name of the solution file
#' @return data.frame with two clumns, idx (index) and group
read_mdgp_solutionfile <- function(file_name) {
  #Sanity check
  stopifnot(file.exists(file_name))

  # Read first line containing n=N and G (no. groups)
  line1 <- readLines(file_name, n=1)
  # Extract number of groups
  n_g <- str_replace(line1, ".*(G = )([0-9]+).*", "\\2") %>% as.numeric()
  # Extract number of individuals n
  n <- str_replace(line1, "^(N = )([0-9]+).*", "\\2") %>% as.numeric()

  #Read the groups
  groups <- read_delim(file=file_name, skip=1 + n_g, delim=" ", col_names=c("id_int", "group")) %>%
    mutate(idx = as.numeric(id_int)+1,
           group = as.numeric(group) + 1) %>%
    select(idx, group)

  return(groups)
}

#' Simple partitioning individuals into groups using simple random sampling
#'
#' @param frame tibble containing a column `id` of the sampling units (e.g. people)
#' @param m requested group size (integer). If m not a divisor of nrow(people) then put leftovers in existing groups
#'
#' @return A list of vectors where each vector containing the ids (row number in frame) of those belonging to the same group

make_partition_srs <- function(frame, m) {
  # Number of groups needed
  n_g <- nrow(frame) %/% m
  # Group membership indicator
  g <- rep(seq_len(n_g), length.out=nrow(frame))
  # Permutation order
  perm <- sample(seq_len(nrow(frame)))

  # Make a list of ids in each group
  groups <- map(seq_len(n_g), function(i) {
    idx <- perm[which(g == i)]
    frame %>% slice(idx) %>% pull(id)
  })

  return(groups)
}


#' Convert a list of sessions into a MDGP distance metric
#'
#' @param past_sessions List of previous sessions
#' @param current_frame The current frame
#'
#' @description A session is a data.frame containing the columns id1, id2 and date
#' @examples
#'  ppl <- data.frame(id=sprintf("id%.3d", 1:100))
#'  session_dates <- seq(as.Date("2021-04-01"), length.out=4, by="1 week")
#'  frames <- map_df( session_dates, ~ ppl %>% slice_sample(n = rbinom(1,nrow(ppl), prob=0.7)) %>% mutate(date=.x))
#'  past_sessions <- frames %>% group_split(date) %>% map(~ rsocialroulette(.x, past_sessions=NULL, m=4, algorithm="srs")) %>% setNames(session_dates) %>% .[1:3]
#'  current_frame <- frames %>% filter(date == max(date))
#'  sessions_to_distance( current_frame, past_sessions) %>% pull(dist) %>% table()
sessions_to_distance <- function(current_frame, past_sessions) {
  # Make all pairs within a group
  group_to_pairs <- function(group) {
    expand_grid(id1=group, id2=group) %>% filter(id1 < id2)
  }

  #Make all potential pairs in current_frame
  current_pairs <- expand_grid(id1=current_frame %>%  pull(id),
                               id2=current_frame %>%  pull(id)) %>%
    filter(id1 < id2) %>%
    mutate(date = current_frame$date[1])

  #If there are no past sessions then no need to do more.
  if (is.null(past_sessions)) {
    return(current_pairs %>% mutate(dist=1))
  }

  #Make the pairs data.frame for each session
  past_pairs <- map_dfr(past_sessions,  ~ map_df(.x, ~ group_to_pairs(.x)), .id="date")

  #Today distance
  Delta <- past_pairs %>% pull(date) %>% unique() %>% as.Date() %>% diff() %>% mean() %>% as.numeric()
  date_no_meet <- (past_pairs %>% pull(date) %>% as.Date() %>%  min() ) - Delta
  dist_today <- ((current_pairs %>% pull(date) %>% .[[1]]) - date_no_meet) %>% as.numeric()

  #Compute dist in days for pairs to past time where they met in a session
  current_dist <- current_pairs %>%
    left_join(past_pairs %>% select(id1, id2, date), by=c("id1", "id2"), suffix=c(".current",".past")) %>%
    mutate(dist = difftime(date.current, date.past, unit="days") %>% round() %>%  as.numeric(),
           dist = if_else(is.na(dist), dist_today, dist))

  #Return
  return(current_dist %>% rename(date = `date.current`) %>% select(-date.past))
}

#' Take a frame including a grouping column and convert this to a partitionList
#'
#' @param frame The frame with `id` and `group` columns to convert
#' @return A list of vectors, where each vector contains all individuals in the corresponding group
frame_to_partitionList <- function(frame) {
  n_g <- length(unique(frame$group))
  map(seq_len(n_g), ~ frame %>% filter(group == .x) %>% pull(id))
}

#' Take a partitionList and convert it to a frame including a grouping column
#'
#' @param A partitionList, i.e. list of vectors, where each vector contains all individuals in the corresponding group
#' @return frame The frame with `id` and `group` columns to convert

partitionList_to_frame <- function(l) {
  map_df( seq_len(length(l)), ~ tibble(id=l[[.x]], group=.x))
}


#' Match a partition based on index (i.e. MGDP output) to ids in the frame
#'
#' @param partition A tibble containing idx and group columns (MDGP output)
#' @param frame A frame (i.e. a tibble with at least a column denoted `id`) to match the idx to. This has to be the frame used to generate the partition in the first place.
#' @return The `frame` tibble augmented with an additional `group` column
partition_to_frame <- function(partition, frame ) {
  #Sanity checks
  stopifnot(nrow(partition) == nrow(frame))

  #Match idx to frame
  frame %>% mutate(idx = row_number()) %>%
    left_join(partition, by="idx") %>%
    select(-idx)
}

#' Make a new lunch roulette session maximizing the gossip to exchange
#'
#' @param current_frame A tibble containing the participants of the current round, i.e. it has a column `id` containing a unique identifier and a `date` column representing the date of the session.
#' @param past_sessions A list of partitionLists, i.e. each partition is a list of vectors containing the id of the members of the corresponding group.
#' @param m minimum group size, i.e. all groups will be at least size m.
#' @return A partitionList containing the partitioning of current_frame maximizing the overall sum of gossip to be exchanged.
#' @export
#' @examples
#' library(socialroulette)
#' #Distribute 5 ppl into groups of at least 2
#' frame <- tibble(id=sprintf("id%.3d@math.su.se", 1:5), date=Sys.Date())
#' rsocialroulette(frame, past_sessions=NULL, m=2, algorithm="srs")
#'
#' #Same, but to test that algorithm="mdgp" also works with past_session=NULL
#' rsocialroulette(frame, past_sessions=NULL, m=2, algorithm="mdgp")
#'
#'# Create a history of previous sessions as well as the current frame from a population of 100
#' ppl <- tibble(id=sprintf("id%.3d@math.su.se", 1:100))
#' session_dates <- seq(as.Date("2021-04-01"), length.out=4, by="1 week")
#' # Simulate changing participation each week (with prob 0.7 they will attend)
#' frames <- map_df( session_dates, ~ ppl %>% slice_sample(n = rbinom(1,nrow(.), prob=0.7)) %>% mutate(date=.x))
#' # Frame for the current session is the last one
#' current_frame <- frames %>% filter(date == max(date))
#' # Simulate the partitions for the last 3 weeks
#' past_sessions <- frames %>% filter(date < max(session_dates)) %>%
#'   group_split(date) %>%
#'   map(~rsocialroulette(current_frame=.x, past_sessions=NULL, m=4)) %>%
#'   setNames(session_dates[1:3])
#' # Make the partition, which tries to bring together those who have not met in a long time
#' rsocialroulette(current_frame, past_session, m=4, algorithm="mdgp")
rsocialroulette <- function(current_frame, past_sessions, m, algorithm=c("mdgp", "srs")) {
  # Sanity checks
  stopifnot( "id" %in% colnames(current_frame))
  stopifnot( "date" %in% colnames(current_frame))
  stopifnot( all(!is.na(as.Date(names(past_sessions)))))
  algorithm <- match.arg(algorithm, choices=c("mdgp", "srs"))

  #Debug info
  cat(str_c("Partitioning ", nrow(current_frame), " individuals into groups of at least ", m, ".\n"))

  #Read output and convert it to a partitionList
  if (algorithm == "mdgp") {
    #Make a specification file and solve it
    spec_file <- write_mdgp_specfile(current_frame, past_sessions, m=m)
    res <- mdgp_solver(spec_file)
    #Read output as a partitionList
    partitionList <- read_mdgp_solutionfile(res$solution_file) %>%
      partition_to_frame(frame=current_frame) %>%
      frame_to_partitionList()
  } else {
    partitionList <- make_partition_srs(current_frame, m=m)
  }

  #Group sizes etc.
  groups <- map_dbl(partitionList, length)

  cat(str_c("Created ", length(groups), " groups of sizes ", str_c(groups, collapse=" "), ".\n"))

  return(partitionList)
}

testit <- function() {
   #library(socialroulette)
   ppl <- tibble(id=sprintf("id%.3d@math.su.se", 1:100))
   session_dates <- seq(as.Date("2021-04-01"), length.out=4, by="1 week")
   frames <- map_df( session_dates, ~ ppl %>% slice_sample(n = rbinom(1,nrow(ppl), prob=0.7)) %>% mutate(date=.x))
   past_sessions <- frames %>% group_split(date) %>% map(~rsocialroulette(.x, past_sessions=NULL, m=4, algorithm="srs")) %>% setNames(session_dates) %>% .[1:3]
   current_frame <- frames %>% filter(date == max(date))

   spec_file <- write_mdgp_specfile(current_frame, past_sessions, m=4)
   res <- mdgp_solver(spec_file)
   current_frame_with_group <- read_mdgp_solutionfile(res$solution_file) %>%
     partition_to_frame(frame=current_frame)

   # Show group sizes
   current_frame_with_group %>% group_by(group) %>% summarise(n=n()) %>%
     pull(n) %>% table()

   #Make a partitionList for the current session
   current_partitionList <-current_frame_with_group %>% frame_to_partitionList()
   current_partitionList

   #Extend the past sessions
   sessions <- past_sessions
   sessions[[current_frame_with_group %>% slice(1) %>% pull(date) %>% as.character()]] <- current_partitionList
   sessions

   #Check if anyone has a distance of zero?
   sessions_to_distance(current_frame, sessions) %>% pull(dist)

   #R CMD INSTALL .
   #Rscript -e "library(socialroulette) ; mdgp()"
}

