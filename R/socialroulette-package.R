#' socialroulette: A package to generate social groupings
#'
#' A package for partitioning individuals into groups of a pre-specified size. This can be as simple as using simple random sampling (srs)
#' to divide n individuals into groups of size at least m. If one keeps track of the past partitions then an additional
#' aim can be to try to maximize the time that people have a reunion, i.e. end up in the same group. This boils
#' down to an instance of the maximally diverse grouping problem.
#' See the [package website](https://hoehleatsu.github.io/socialroulette/)  for more information, documentation and examples.
#'
#' A side effect of the package is that it provides a maximally diverse grouping problem solver, which
#' can in principle be used for other purposes. As a consequence, internal functionality is exposed using
#' export statements.
#'
#' @docType package
#' @keywords package
#' @name socialroulette-package
#' @importFrom stats dist
#' @importFrom stats setNames
#' @importFrom tibble tibble
#' @importFrom readr read_delim write_delim
#' @import dplyr
#' @import purrr
NULL

