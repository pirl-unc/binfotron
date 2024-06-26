% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stackable_functions.R
\name{round_data}
\alias{round_data}
\title{Round data.table columns}
\usage{
round_data(
  my_dt = NULL,
  col_names = NULL,
  digits = 0,
  my_summary = "Rounding data.",
  sample_key = get_default_sample_key(),
  readme_path = NULL
)
}
\arguments{
\item{my_dt}{data.table input}

\item{col_names}{Vector of character strings to name the columns that will have this 
operation performed on them.  Uses \code{\link{operatable_columns}}}

\item{digits}{Integer to specify how mnay decimal digits to round to.
See \code{\link{round}} 'digits' argument.}

\item{my_summary}{Character string to change the default comment that will be appended 
to my_dt.}

\item{sample_key}{Character string to specify the column that is the sample key. This 
column will not be operated upon.}

\item{readme_path}{Optional path to which the comments will be appended.}
}
\description{
A stackable wrapper for \code{\link{round}} to turn it into a stackable transformation
}
\details{
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
round_data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}
