

#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' import_model_data
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Perform common steps needed for importing model data
#' 
#' @description
#' All models will start with a matrix of ranked features, and the feature
#' values of the features. This fxn orders the features, pulls out the top x and
#' outputs disc and validation dataset matrices with the correct indep and dep
#' variables. This function has the option to specify certain features or to pull
#'the top x feature from a columns ranks file.
#'  
#' @param import_paths Vector of paths to import and merge
#' @param ranks_feature_path Path to data matrix with feature names in one column and ranks in another.
#' @param ranks_clm String to indicate which column contains the column ranks
#' @param ranks_count Integer to indicate how many features to take from the column ranks
#' @param merge_file_patterns Character vector of file patterns to look for and merge to make the final output matrices
#' @param set_clm Column name to indicate which samples are in "Validation" vs "Discovery" sets
#' @param sample_key Column name to use as the sample key to merge all matrices
#' @param indep_vars Character string of column names for the independent variables
#' @param dep_vars Character string of column names for the dependent variables
#' 
#' @export
import_model_data = function(
	import_paths,
	ranks_feature_path = NULL,
	ranks_clm = NULL,
	ranks_count = NULL,
	set_clm = "Set",
	sample_key = binfotron::get_default_sample_key(),
	indep_vars = NULL,
	dep_vars = NULL
){
	
	if (!is.null(indep_vars) & (!is.null(ranks_feature_path) ))
		stop("Ranks arguments and indep_vars should not be set at the same time.")
	
	if (!is.null(ranks_feature_path)){
		indep_vars = get_features_from_ranks(ranks_feature_path, ranks_count, ranks_clm)
	}
	
	cat(paste0("dep_vars: ", paste0( dep_vars, collapse= ", "), "\n"))
	cat(paste0("indep_vars: ", paste0( indep_vars, collapse= ", "), "\n"))
	import_clms = c(sample_key, set_clm, dep_vars, indep_vars)
	if (any(duplicated(import_clms))){
		duplicated_clms = import_clms[duplicated(import_clms)]
		message(paste0("duplicated_clms: ", paste0( duplicated_clms, collapse= ", ")))
		stop("You have duplicated column names specified in your set_clm, sample_key, indep_vars, and dep_vars. Please resolve this.")
	}
	
	# limit to variables that are in import_clms
	import_df = serial_import(
		import_paths = import_paths,
		sample_key = sample_key,
		import_clms = import_clms
	)
		
		
	if (!is.null(set_clm)){
		if (!all(c("Validation","Discovery") %in% import_df[[set_clm]])) stop(paste0("Your set_clm, ", set_clm, ", does not contain 'Discovery' and 'Validation' strings."))
		val_df = import_df[import_df[[set_clm]] == "Validation", ]
		import_df = import_df[import_df[[set_clm]] == "Discovery", ]
		val_df = val_df[, names(val_df)[names(val_df) != set_clm]]
		import_df = import_df[, names(import_df)[names(import_df) != set_clm]]
	} else {
		val_df = NULL
	}
	
	# dataset prep is the right place to make sure the features aren't univariable for discovery validation sets, not here
	return(list(disc_df=import_df, val_df=val_df, indep_vars=indep_vars))
}


#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' get_features_from_ranks
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Return vector of features names from a data.frame with rank values.
#'  
#' @param ranks_feature_path Path to data matrix with feature names in one column and ranks in another.
#' @param ranks_clm String to indicate which column contains the column ranks. \code{NULL} will return the second column.
#' @param ranks_count Integer to indicate how many features to take from the column ranks, \code{NULL} will return all the features ordered by rank.
#' 
#' @export
get_features_from_ranks = function(
	ranks_feature_path,
	ranks_count = NULL,
	ranks_clm = NULL
){
	ranks_df = fread(ranks_feature_path, data.table = F)
	if (is.null(ranks_clm)) ranks_clm = names(ranks_df)[2]
	ranks_cls = class(ranks_df[[ranks_clm]])
	if (ranks_cls %in% c('numeric', 'integer')){
		ranks_df = ranks_df[order(ranks_df[[ranks_clm]]),]
	} else if (ranks_cls %in% c('logical')) {
		ranks_df = ranks_df[rev(order(ranks_df[[ranks_clm]])),]
	} else {
		stop(paste0("Unrecognized class for the ", ranks_clm," ranks_clm: ", ranks_cls))
	}
	if (is.null(ranks_count)){
		indep_vars = ranks_df[[1]]
	} else {
		indep_vars = ranks_df[[1]][1:ranks_count]
	}
	return(indep_vars)
}


#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' serial_import
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Import a vector fo file paths into one matrix.
#'  
#' @param import_paths Vector of paths to import and merge
#' @param sample_key Column name to use as the sample key to merge all matrices
#' @param import_clms Vector of column names to import
#' 
#' @export
serial_import = function(
  import_paths,
  sample_key = binfotron::get_default_sample_key(),
  import_clms = NULL
){
	import_df = NULL
	for (import_path in import_paths){
		this_df = fread(import_path, data.table = F)
		if (!sample_key %in% names(this_df)) stop(paste0("This path does not contain your sample_key, ", sample_key,": ", import_path))
		if(!is.null(import_clms)) this_df = this_df[, import_clms[import_clms %in% names(this_df)]]
		if (ncol(this_df) == 1) stop(paste0("This file did not have any data you were looking for: ", import_path))
		if(is.null(import_df)){
			import_df = this_df
		} else {
			import_df = merge(import_df, this_df, by = sample_key, all = T)
		}
	}
	
	if(!is.null(import_clms)){
		missing_clms = import_clms[!(import_clms %in% names(import_df))]
		if (length(missing_clms) > 0) stop(paste0(c("Your imported data are missing the following columns:", missing_clms), collapse= "\n"))
	}
	return(import_df)
}
