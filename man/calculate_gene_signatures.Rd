% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stackable_functions.R
\name{calculate_gene_signatures}
\alias{calculate_gene_signatures}
\title{Calculate gene signatures on my_dt using a gmt file.}
\usage{
calculate_gene_signatures(
  my_dt = NULL,
  gene_cols = NULL,
  gene_element = 1,
  gmt_file_path,
  min_genes = 1,
  my_fun = mean,
  my_summary =
    "Gene signatures by default are the mean of the genes listed in the signature for each sample.",
  only_return_signatures = TRUE,
  summary_output_path = NULL,
  sample_key = get_default_sample_key(),
  signatures = NULL,
  readme_path = NULL,
  na_remove = TRUE
)
}
\arguments{
\item{my_dt}{data.table input}

\item{gene_cols}{Character vector indicating the columns that should be used in the 
gene signature. Uses \code{operatable_columns}.  If \code{NULL} if will take all 
numeric columns.}

\item{gene_element}{For gene names that have a pipe in them, which position should be used ("1|2 etc"). Integer.}

\item{gmt_file_path}{String path to the gmt file.}

\item{min_genes}{Integer for the minimum number of genes that have to be found from a
gene signature to keep it.}

\item{my_fun}{Function used to combine the genes for each sample. Default is \code{\link{mean}}}

\item{my_summary}{Character string to change the default comment that will be appended 
to my_dt.}

\item{only_return_signatures}{Boolean indicating whether or not all of the non-gene columns
should be included in the result.  Gene columns will be dropped.  This is helpful in 
returning categories with the gene signature scores}

\item{summary_output_path}{If specified this is where the summary data for making this gene 
signature will go.}

\item{sample_key}{Character string to specify the column that is the sample key. This 
column will not be operated upon.}

\item{signatures}{Specifies which gene signtures on the gmt file will be used. If \code{NULL}
all of the signatures wil be used.}

\item{readme_path}{Optional path to which the comments will be appended.}

\item{col_names}{Vector of character strings to name the columns that will have this 
operation performed on them.  Uses \code{\link{operatable_columns}}}
}
\description{
Calculates gene signatures on my_dt using a gmt file. The method of calculating the 
gene signture can be changed using \code{my_fun}
}
\details{
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
calculate_gene_signatures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}
\section{Todos}{

\itemize{
  \item Should convert a data.frame to a data.table and then convert it back for the output.
}
}

\seealso{
Other gene_signature: 
\code{\link{bgv_lab_prep_count_data_for_gene_signature}()}
}
\concept{gene_signature}
