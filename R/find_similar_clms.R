#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' find_similar_clms
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' @title Finds columns that have a high degree of correlation with this_clm
#' 
#' @description
#' Takes \code{this_clms} and compares it to all other columns of \code{my_df} using \code{\link{stats::cor.test}}.  If the absolute value of that row is greater than \code{acceptable_similarity}, then the either a boolean is returned true or the similar columns are returned.
#' 
#' @param this_clm string to specify column that will be compared to the others
#' @param my_df \code{data.frame} to search
#' @param test_clms character vector of columns to test. If left blank all, but \code{this_clm} and \code{my_key} will be checked.
#' @param acceptable_similarity number to indicate the max absolute rho value that will be allowed with other columns.
#' @param return_boolean boolean to specify if a boolean is returned or the \code{found_clms}
#' @param my_key string to specify key column of the \code{data.frame}
#' @param corr_method  string to specify the test used for the correlation. Passed to \code{method} arg of \code{\link{stats::cor.test}}
#' 
#' @export
find_similar_clms = function(
	this_clm,
	my_df,
	test_clms = NULL,
	acceptable_similarity = 0.9,
	return_boolean = T, # faster in that it can quit after finding the first problem
	my_key = get_default_sample_key(),
	corr_method = 'spearman'
){
	if(is.null(test_clms)) test_clms = names(my_df)[ ! names(my_df) %in% c( this_clm, my_key ) ]
	found_clms =  character(0)
	for (test_clm in test_clms) {
		sim_df = my_df[,c(test_clm, this_clm)]
		sim_df = sim_df[complete.cases(sim_df),]
		if (nrow(sim_df) > 3)
			tryCatch({
				my_cor = cor.test(x=as.numeric( sim_df[[ test_clm ]] ), y=as.numeric( sim_df[[ this_clm ]] ), method = corr_method, exact = FALSE)
				abs_rho = abs(my_cor$estimate)
				if ( !is.na(abs_rho) & abs_rho > acceptable_similarity ) {
					message( test_clm, " is too similar to ", this_clm)
					found_clms = c( found_clms, test_clm )
					if (return_boolean) break
				}
			}, error=function(error_message) {
				message(error_message)
			})
	}
	
	if (return_boolean){
		return( length(found_clms) > 0 )
	} else {
		return( found_clms )
	}
}