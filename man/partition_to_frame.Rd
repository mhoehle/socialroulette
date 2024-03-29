% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/socialroulette.R
\name{partition_to_frame}
\alias{partition_to_frame}
\title{Take a partition and convert it to frame representation}
\usage{
partition_to_frame(l)
}
\arguments{
\item{l}{A partition, i.e. list of vectors, where each vector contains all individuals in the corresponding group}
}
\value{
frame The frame with `id` and `group` columns to convert
}
\description{
The resulting \code{tibble} will have an additional column \code{group}
}
\examples{
round1 <- list(c("id02", "id03", "id04"), c("id05", "id01"))
socialroulette::partition_to_frame(round1)
}
