#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' find_low_variance_clms
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' @title Returns column names of those whose most common factor is used by more than \code{max_fraction_largest} of the samples.
#' 
#' @description
#' Checks all columns of \code{my_df} except \code{my_key} and checks if the most common factor is more than \code{max_fraction_largest}.
#' 
#' @param my_df \code{data.frame} to search
#' @param max_fraction_largest number to indicate the max proportion the top fraction can take up before being reported.
#' @param my_key string to specify key column of the \code{data.frame}
#' 
#' @export
find_low_variance_clms = function(
	my_df,
	max_fraction_largest = max_fraction_largest,
	my_key = get_default_sample_key()
){
	potential_drop_clms = character(0)
	for(clm_num in 1:ncol(my_df)){
		# clm_num = 2
		clm_name = names(my_df)[clm_num]
		if(clm_num != my_key){
			my_values = trimws(my_df[[clm_name]])
			my_values = my_values[!is.na(my_values)]
			my_values = my_values[my_values != ""]
			max_rep = max(summary(factor(my_values), na.rm = T))/length(my_values)
			if ( max_rep > max_fraction_largest ) {
				a(paste0(clm_name, ": largest factor ",max_rep,"\n"))
				potential_drop_clms = c(potential_drop_clms, clm_name)
			}
		}
	}
	return(unique(potential_drop_clms))
}

