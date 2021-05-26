#' Maximally diverse grouping problem solver based on Lai and Hao (2016)
#'
#' Given a distance metric for each possible pair in the frame of n individuals,
#' find a partition into as many groups as possible with group size at least m. This
#' is an instance of the maximallz diverse grouping problem, which we use the algorithm
#' by Lai and Hao (2016) to solve.
#'
#' The maximally diverse grouping problem is about partitioning \eqn{n} individuals
#' into groups of size at least \eqn{m} while maximizing a sum of utility values computed
#' by summing the utility \eqn{u(i,j)} over all individuals \eqn{i},\eqn{j} partitioned into the same group.
#' and summing this quantity over all groups. More formally, let \eqn{d_{ij}} denote the number of
#' time unit (typically days) ago, that individual \eqn{i} and \eqn{j} were in the same group.
#' Note: This distance is derived by looking at the previous partitions and it is a matter
#' of definition what this value should be, if \eqn{i} and \eqn{j} have not previously been
#' in the same group. Let \eqn{G=n\> \text{div}\> m} denote the resulting number of groups where \eqn{\text{div}} denotes integer division.
#' For a given partition let \eqn{x_{ig}} be an indicator variable, which is 1, if \eqn{i} is assigned into
#' group \eqn{g} and zero otherwise.
#' A solver of the maximally diverse grouping problem now tries to maximize
#' \deqn{\sum_{g=1}^G \sum_{i=1}^n \sum_{j=i+1}^n d_{ij} x_{ig} x_{jg},}
#' subject to the conditions
#' \deqn{\sum_{g=1}^{G} x_{ig} = 1, \quad i=1,\ldots,n,}
#' \deqn{\sum_{i=1}^n x_{ig} = n_g, \quad g=1,\ldots,G,}
#' \deqn{x_{ig} \in \{0,1\}, \quad i=1,\ldots,n; g=1,\ldots,G,}
#' where \eqn{n_g} is the size of group \eqn{g}, i.e. \eqn{\sum_{g=1}^G n_g = n}. We shall adopt the convention
#' that group sizes are determined by assigning the labels \eqn{1,\ldots,G} to all individuals and then count
#' how often each label occurs. This means, e.g., that for \eqn{n=7} and \eqn{m=4} we get \eqn{n_g=7\> \text{div}\> 4=1} group with 7 members.
#'
#' Note: The code calls C++ code by Xiangjing Lai and Jin-Kao Hao,
#' for this several temporary files are generated using `tempfile()`.
#'
#' @param mdgp_format_file Path to the file containing the MDGP specification in mdgplib format
#' @param time_limit Number of seconds to iteratively optimize each run. The larger the number of participants to group, the larger this value should be. Rule of thumb: time_limit = exp(0.5 + 0.0025*n)
#' @return File name of the solution file
#'
#' @references Xiangjing Lai and Jin-Kao Hao (2016). *Iterated maxima search for the maximally
#' diverse grouping problem*. European Journal of Operational Research, 254(3), pp. 780-800,
#' https://doi.org/10.1016/j.ejor.2016.05.018
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

#' Convert partition specification to Lai and Hao (2016) MDGP solver input specification.
#'
#' Generate a \code{tempfile()} containing the appropriate specification formed from \code{current_frame} and \code{past_partitions}.
#'
#' The specification format has the following header line
#' \deqn{n\> m\> ds\> n_1\> n_1\> \ldots n_G \> n_G}
#' Here, each group size \eqn{n_g, g=1,\ldots,G} is repeated twice, because the solver has the opportunity to specify a lower as well as an
#' upper boundary for each group size. The keyword `ds` for the solver in the above means that the groups can be of different sizes as determined by the respective \eqn{n_g}.
#'
#' @param current_frame A pair-distance tibble. The frame needs to contain a column denoted `id`
#' @param past_partitions A named list of partitions, where the names correspond to the date when the partition was used.
#' @param m The minimum group size, i.e. all groups have size at least m
#' @return File name of the generated specification file
#' @examples
#' frame  <- tibble::tibble(id=sprintf("id%.2d", 1:5), date=as.Date("2021-04-28"))
#' past_partitions <- list("2021-04-21"=list(c("id02", "id03", "id04"), c("id01", "id05")))
#' spec_file <- socialroulette:::mdgp_write_specfile(frame, past_partitions=past_partitions, m=2)
#' cat(stringr::str_c(readLines(spec_file), collapse="\n"))
mdgp_write_specfile <- function(current_frame, past_partitions, m) {
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
  dist_frame <- partitions_to_distance(current_frame, past_partitions)

  # Make alternative index number (starting from zero) for MDGP-solver instead of the id
  ids <- current_frame %>% dplyr::mutate(idx = dplyr::row_number() - 1) %>% dplyr::select(id, idx)
  # and add this to the dist_frame
  dist_frame <-  dplyr::left_join(dist_frame, ids, by=c("id1"="id")) %>%
    dplyr::left_join(ids, by=c("id2"="id")) %>%
    dplyr::rename(idx.id1 = dplyr::all_of("idx.x"), idx.id2=dplyr::all_of("idx.y")) %>%
    #Make sure idx.id1 < idx.id2 as this is needed for the solver
    dplyr::mutate( idx.id1.ord = pmin(idx.id1, idx.id2),
                   idx.id2.ord = pmax(idx.id1, idx.id2)) %>%
    dplyr::arrange(all_of("idx.id1.ord"), dplyr::all_of("idx.id2.ord"))

  #Sanity check
  stopifnot( all(dist_frame$idx.id1.ord < dist_frame$idx.id2.ord))

  #Make a temporary file
  tmp_file <- str_c(tempfile(), ".txt")
  writeLines(spec, tmp_file)
  readr::write_delim(dist_frame %>% dplyr::select(idx.id1.ord,idx.id2.ord,dist), file=tmp_file, col_names=FALSE, append=TRUE)

  #Done - return filename
  return(tmp_file)
}

#' Read output from the Lai and Hao (2016) MDGP solver
#' @param file_name Name of the solution file
#' @return data.frame with two clumns, idx (index) and group
mdgp_read_solutionfile <- function(file_name) {
  #Sanity check
  stopifnot(file.exists(file_name))

  # Read first line containing n=N and G (no. groups)
  line1 <- readLines(file_name, n=1)
  # Extract number of groups
  n_g <- str_replace(line1, ".*(G = )([0-9]+).*", "\\2") %>% as.numeric()
  # Extract number of individuals n
  n <- str_replace(line1, "^(N = )([0-9]+).*", "\\2") %>% as.numeric()

  #Read the groups (col_types argument used to avoid output - see ?readr::read_delim)
  groups <- readr::read_delim(file=file_name, skip=1 + n_g, delim=" ", col_names=c("id_int", "group"), col_types = cols()) %>%
    dplyr::mutate(idx = as.numeric(id_int)+1,
           group = as.numeric(group) + 1) %>%
    dplyr::select(idx, group)

  return(groups)
}

#' Simple partitioning individuals into groups using simple random sampling
#'
#' @param frame tibble containing a column `id` of the sampling units (e.g. people)
#' @param m requested group size (integer). If m not a divisor of nrow(people) then put leftovers in existing groups
#' @return A list of vectors where each vector containing the ids (row number in frame) of those belonging to the same group
#'
#' @keywords internal

make_partition_srs <- function(frame, m) {
  # Number of groups needed
  n_g <- nrow(frame) %/% m
  # Group membership indicator
  g <- rep(seq_len(n_g), length.out=nrow(frame))
  # Permutation order
  perm <- sample(seq_len(nrow(frame)))

  # Make a list of ids in each group
  groups <- purrr::map(seq_len(n_g), function(i) {
    idx <- perm[which(g == i)]
    frame %>% dplyr::slice(idx) %>% dplyr::pull(id) %>% sort()
  })

  return(groups)
}

#' Make all pairs within a group
#'
#' @param group A vector of ids belonging to the same group
#' @return a data.frame containing all pairs with id of the first being smaller than the id of the second entry
#' @keywords internal
group_to_pairs <- function(group) {
  tidyr::expand_grid(id1=group, id2=group) %>% dplyr::filter(id1 < id2)
}

#' Convert a list of partitions into a data.frame with all pairs
#'
#' The name of the list contain the date at which each partition was used.
#' @param partitions A list of past partitions, i.e. a named list where each entry is a partition and the names reflect the dates the partition were used
#' @export
#' @examples
#' partitions <- list("2021-04-21"=list(c("id02", "id03", "id04"), c("id05", "id01")),
#'                   "2021-04-28"=list(c("id05", "id03", "id04"), c("id02", "id01")))
#' socialroulette::partitions_to_pairs(partitions)

partitions_to_pairs <- function(partitions) {
  purrr::map_dfr(partitions,  ~ purrr::map_df(.x, ~ group_to_pairs(.x)), .id="date")
}

#' Function to convert a data.frame of pairs to a partition
#'
#' @param pairs_df data.frame containing the pairs, i.e. it has columns date, id1 and id2
#' @return A named list of partitions
#' @export
#' @keywords internal
pairs_to_partition <- function(pairs_df) {
  res <- list()
  for (i in seq_len(nrow(pairs_df))) {
    ids <- pairs_df[i, c("id1", "id2")] %>% as.character()
    in_list <- purrr::map_lgl(res, function(l) any(ids %in% l))
    if (any(in_list)) {
      if (sum(in_list) == 1) {
        # Add to bucket
        res[[which(in_list)]] <- c(res[[which(in_list)]], ids) %>% unique() %>% sort()
      } else {
        # Merge two buckets
        res[[length(res)+1]] <- res[which(in_list)] %>% unlist %>% unique() %>% sort()
        res[which(in_list)] <- NULL
      }
    } else {
      # Make a new bucket
      res[[length(res)+1]] <- ids
    }
  }
  return(res)
}

#' Convert a pairs data.frame into a partition list
#'
#' The function works by splitting the pairs data.frame by date and then apply the internal
#' function pairs_to_partition to the resulting sub data.frame.
#'
#' @param pairs A \code{data.frame} with columns \code{date} and \code{id1} and \code{id2}.
#' @seealso pairs_to_partition
#' @export
#' @examples
#' partitions <- list("2021-04-21"=list(c("id02", "id03", "id04"), c("id01", "id05")),
#'                   "2021-04-28"=list(c("id03", "id04", "id05"), c("id01", "id02")))
#' partitions2 <- partitions %>% socialroulette::partitions_to_pairs() %>%
#'                           socialroulette::pairs_to_partitions()
#' all.equal(partitions, partitions2)
pairs_to_partitions <- function(pairs) {

  partitions <- pairs %>% dplyr::arrange(date) %>% dplyr::group_split(date) %>%
    purrr::map(~ pairs_to_partition(.x)) %>%
    setNames(pairs %>% dplyr::pull(date) %>% unique())

  return(partitions)
}


#' Convert a list of partitions into pairs format including a distance metric
#'
#' @param past_partitions List of previous partitions
#' @param current_frame The current frame
#' @export
#' @examples
#' frame  <- tibble::tibble(id=sprintf("id%.2d", 1:5), date=as.Date("2021-04-28"))
#' partitions <- list("2021-04-14"=list(c("id02", "id03", "id04"), c("id01", "id05")),
#'                   "2021-04-21"=list(c("id03", "id04", "id05"), c("id01", "id02")))
#' socialroulette::partitions_to_distance(frame, partitions)

partitions_to_distance <- function(current_frame, past_partitions) {
  #Make all potential pairs in current_frame
  current_pairs <- tidyr::expand_grid(id1=current_frame %>%  dplyr::pull(id),
                               id2=current_frame %>%  dplyr::pull(id)) %>%
    dplyr::filter(id1 < id2) %>%
    dplyr::mutate(date = current_frame$date[1])

  #If there are no past partitions then no need to do more.
  if (is.null(past_partitions)) {
    return(current_pairs %>% dplyr::mutate(dist=1))
  }

  #Make the pairs data.frame for each partition
  past_pairs <- partitions_to_pairs( past_partitions)

  #Define distance to today for those who have not met so far - special case if only one past value
  Delta <- rbind(current_pairs, past_pairs) %>%
    dplyr::pull(date) %>% unique() %>%
    as.Date() %>%
    sort() %>%
    diff() %>%
    mean() %>% as.numeric()
  # Compute with Delta
  date_no_meet <- (past_pairs %>% dplyr::pull(date) %>% as.Date() %>%  min() ) - Delta
  dist_today <- ((current_pairs %>% dplyr::pull(date) %>% .[[1]]) - date_no_meet) %>% as.numeric()

  #Compute dist in days for pairs to past time where they met in a partition
  current_dist <- current_pairs %>%
    dplyr::left_join(past_pairs %>% dplyr::select(id1, id2, date), by=c("id1", "id2"), suffix=c(".current",".past")) %>%
    dplyr::mutate(dist = difftime(date.current, date.past, units="days") %>% round() %>%  as.numeric(),
           dist = dplyr::if_else(is.na(dist), dist_today, dist))

  #Return
  return(current_dist %>% dplyr::rename(date = `date.current`) %>% dplyr::select(-date.past))
}

#' Take a frame and convert this to a corresponding partition
#'
#' This is a small useful converter function for taking a tibble including a `group` and `id` column and convert it to a partition.
#'
#' @param frame The frame with `id` and `group` columns to convert
#' @return A partition, i.e. a list of vectors, where each vector contains all individuals in the corresponding group
#' @keywords internal
#' @export
#' @examples
#' p <- tibble::tibble(id=sprintf("id%.02d",1:5), group=c(1,1,1,2,2))
#' socialroulette::frame_to_partition(p)
frame_to_partition <- function(frame) {
  n_g <- length(unique(frame$group))
  purrr::map(seq_len(n_g), ~ frame %>% filter(group == .x) %>% dplyr::pull(id) %>% sort())
}

#' Take a partition and convert it to frame representation
#'
#' The resulting \code{tibble} will have an additional column \code{group}
#'
#' @param l A partition, i.e. list of vectors, where each vector contains all individuals in the corresponding group
#' @return frame The frame with `id` and `group` columns to convert
#' @export
#' @examples
#' round1 <- list(c("id02", "id03", "id04"), c("id05", "id01"))
#' socialroulette::partition_to_frame(round1)
partition_to_frame <- function(l) {
  purrr::map_df( seq_len(length(l)), ~ tibble::tibble(id=l[[.x]], group=.x))
}


#' Match a partition based on index (i.e. MGDP output) to ids in the frame
#'
#' @param mdgp_partition A tibble containing idx and group columns (MDGP output)
#' @param frame A frame (i.e. a tibble with at least a column denoted `id`) to match the idx to. This has to be the frame used to generate the partition in the first place.
#' @return The `frame` tibble augmented with an additional `group` column
#' @keywords internal
mdgp_partition_to_frame <- function(mdgp_partition, frame ) {
  #Sanity checks
  stopifnot(nrow(mdgp_partition) == nrow(frame))

  #Match idx to frame
  frame %>% dplyr::mutate(idx = dplyr::row_number()) %>%
    dplyr::left_join(mdgp_partition, by="idx") %>%
    dplyr::select(-idx)
}

#' Make a new lunch roulette partition maximizing the gossip to exchange
#'
#' One can either use simple random sampling (srs) to generate the partition or
#' solve the maximally diverse grouping problem (mdgp) using the algorithm by
#' Lai and Hao (2016) in order to maximize time since last meets over all groups.
#'
#' @param current_frame A tibble containing the participants of the current round, i.e. it has a column `id` containing a unique identifier and a `date` column representing the date of the partition.
#' @param past_partitions A list of partition lists each named by the date the partition was used (this is used to determine temporal distance to last meeting)
#' @param m minimum group size, i.e. all groups will be at least size m.
#' @param algorithm String specifying either "srs" (simple random sampling - DEFAULT) or "mdgp" (maximally diverse grouping problem)
#' @param \dots Additional arguments to be sent to the solver
#' @return A partitioning of current_frame maximizing the overall sum of gossip to be exchanged.
#'
#' @seealso mdgp_solver
#' @references Höhle M (2021), Long time, no see: Virtual Lunch Roulette, Blog post, \url{https://staff.math.su.se/hoehle/blog/2021/04/04/socialsamp.html}
#' @references Xiangjing Lai and Jin-Kao Hao (2016). *Iterated maxima search for the maximally
#' diverse grouping problem*. European Journal of Operational Research, 254(3), pp. 780-800,
#' https://doi.org/10.1016/j.ejor.2016.05.018
#' @examples
#' today <- Sys.Date()
#' frame <- tibble::tibble( id=sprintf("id%.02d",1:5), date=today)
#' round1 <- rsocialroulette(current_frame = frame, m=2, algorithm="srs")
#' round1
#'
#' #Generate list of past partitions
#' past_partitions <- list(round1) %>% setNames(today)
#' frame2 <- frame %>% dplyr::mutate(date = today+7)
#' round2 <- rsocialroulette(current_frame = frame2,
#'                           past_partitions=past_partitions, m=2, algorithm="mdgp")
#' round2
#' @export
rsocialroulette <- function(current_frame, past_partitions=NULL, m, algorithm=c("mdgp", "srs"), ...) {
  # Sanity checks
  stopifnot( "id" %in% colnames(current_frame))
  stopifnot( "date" %in% colnames(current_frame))
  stopifnot( all(!is.na(as.Date(names(past_partitions)))))
  algorithm <- match.arg(algorithm, choices=c("mdgp", "srs"))
  if (!is.null(past_partitions) & (algorithm == "srs")) {
    stop("Simple random sampling (srs) does not work with past partition information.")
  }
  #Debug info
  cat(str_c("Partitioning ", nrow(current_frame), " individuals into groups of at least ", m, " (", ifelse(is.null(past_partitions),"no past partitions", "Optimizing wrt. past partitions"),").\n"))

  #Read output and convert it to a partition
  if (algorithm == "mdgp") {
    #Make a specification file and solve it
    spec_file <- mdgp_write_specfile(current_frame, past_partitions, m=m)
    res <- mdgp_solver(spec_file, ...)
    #Read output as a partition
    partition <- mdgp_read_solutionfile(res$solution_file) %>%
      mdgp_partition_to_frame(frame=current_frame) %>%
      frame_to_partition()
  } else {
    partition <- make_partition_srs(current_frame, m=m)
  }

  #Group sizes etc.
  groups <- purrr::map_dbl(partition, length)

  cat(str_c("Created ", length(groups), " groups of sizes ", str_c(groups, collapse=" "), ".\n"))

  return(partition)
}
