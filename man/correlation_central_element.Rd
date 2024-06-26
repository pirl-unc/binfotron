% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_central_clone.R
\name{correlation_central_element}
\alias{correlation_central_element}
\title{Determine the most similar of a list of elements to all other elements in the list based on a correlation matrix}
\usage{
correlation_central_element(correlation_data, elements)
}
\arguments{
\item{correlation_data}{Data.frame or matrix of similarity scores / correlations containing at a minimum all elements. This df should be symmetrical - i.e. have equivalent values for row A to column B and row B to column A}

\item{elements}{Character vector of the set of elements ( i.e. features ) from which to determine greatest similarity. Values must be a subset of correlation_data row and column names}
}
\value{
Returns the name of the element that is most similar to all other elements in the set
}
\description{
Calculates maximum similarity/correlation between any single element of group and all other elements
}
