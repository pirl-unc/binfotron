% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_central_clone.R
\name{write_ranked_central_elements_table}
\alias{write_ranked_central_elements_table}
\title{Saves to file a dataframe with two columns, Element and Rank}
\usage{
write_ranked_central_elements_table(
  central_elements,
  rank_data,
  output_path = file.path(getwd(), "ranked_unique_central_elements_ranks.tsv"),
  display_output = F
)
}
\arguments{
\item{central_elements}{Named list of central element from each cluster ( as returned by isolate_central_cluster_elements ), assumes naming convention of {cluster_group}_{cluster_number_within_group}}

\item{output_path}{One-length character vector with fully qualified path to output file}

\item{display_output}{Boolean option to display final output}

\item{element_ranks}{Named integer vector of element rankings}
}
\value{
Returns data.frame of unique ranked central cluster elements as saved to file
}
\description{
Creates a file with one tab separated row per each unique central cluster element along with its ranking
}
