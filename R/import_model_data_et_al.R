#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' import_model_data
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Perform common steps needed for importing model data.
#' 
#' @description
#' All models will start with a matrix of ranked features, and the feature
#' values of the features. This fxn orders the features, pulls out the top x and
#' outputs disc and validation dataset matrices with the correct indep and dep
#' variables. If you don't need to get features from ranks, import paths and
#' split into discovery datasets just use the individual pieces:  
#' get_top_ranked_features, serial_merge, split_disc_val
#'  
#' @param ranks_feature_path Path to file that contains the ranking of features
#' @param ranks_clm String to indicate which column contains the column ranks
#' @param ranks_count Integer to indicate houw many features to take from the column ranks
#' @param merge_paths Character vector of file paths to merge into one matrix
#' @param set_clm Column name to indicate which samples are in "Validation" vs "Discovery" sets
#' @param sample_key Column nameto use as the sample key to merge all matrices
#' @param indep_vars Character string of column names for the independent variables
#' @param dep_vars Character string of column names for the dependent variables
#' 
#' @export
import_model_data = function(
	merge_paths,
	ranks_feature_path = NULL,
	ranks_clm = "Ranks",
	ranks_count = NULL,
	set_clm = "Set",
	sample_key = binfotron::get_default_sample_key(),
	indep_vars = NULL,
	dep_vars = NULL
){
	
	if (!is.null(indep_vars) & (!is.null(ranks_feature_path) | !is.null(ranks_count) | !is.null(ranks_clm)))
		stop("Ranks arguments and indep_vars should not be set at the same time.")
	
	if (!is.null(ranks_count)){
		if (is.null(ranks_feature_path)){
			stop("ranks_feature_path cannot be NULL if using ranks_count.")
		} else {
			ranks_df = fread(ranks_feature_path, data.table = F)
			indep_vars = get_top_ranked_features(ranks_feature_path, ranks_count, ranks_clm)
		}
	}
	
	cat(paste0("dep_vars: ", paste0( dep_vars, collapse= ", "), "\n"))
	cat(paste0("indep_vars: ", paste0( indep_vars, collapse= ", "), "\n"))
	all_clms = c(sample_key, set_clm, dep_vars, indep_vars)
	if (any(duplicated(all_clms))){
		duplicated_clms = all_clms[duplicated(all_clms)]
		message(paste0("duplicated_clms: ", paste0( duplicated_clms, collapse= ", ")))
		stop("You have duplicated column names specified in your set_clm, sample_key, indep_vars, and dep_vars. Please resolve this.")
	}
	
	# limit to variables that are in all_clms
	import_df = serial_merge(merge_paths, sample_key, all_clms)

	if (!is.null(set_clm)){
		split_result = split_disc_val(import_df,set_clm)
		import_df = split_result$disc_df
		val_df = split_result$val_df
		rm(split_result)
	} else {
		val_df = NULL
	}
	
	# dataset prep is the right place to make sure the features aren't univariable for discovery validation sets, not here
	return(list(disc_df=import_df, val_df=val_df, indep_vars=indep_vars))
}




#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' get_top_ranked_features
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Output the top ranked features 
#'
#' @param ranks_df data.frame that includes the feature names and column by which to rank them
#' @param ranks_clm String to indicate which column contains the column ranks
#' @param ranks_count Integer to indicate how many features to take from the column ranks
#' @export
function get_top_ranked_features(
	ranks_df,
	ranks_count,
	ranks_clm = "Ranks",
){
	ranks_cls = class(ranks_df[[ranks_clm]])
	if (ranks_cls %in% c('numeric', 'integer')){
		ranks_df = ranks_df[order(ranks_df[[ranks_clm]]),]
	} else if (ranks_cls %in% c('logical')) {
		ranks_df = ranks_df[rev(order(ranks_df[[ranks_clm]])),]
	} else {
		stop(paste0("Unrecognized class for the ", ranks_clm," ranks_clm: ", ranks_cls))
	}
	top_features = ranks_df[[1]][1:ranks_count]
	return(top_features)
}

#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' serial_merge
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Merge all matrices from a vector of paths
#'
#' @param merge_paths Vector of paths to data matrices
#' @param sample_key String to indicate which column contains the sample keys
#' @param import_clms Vector of column names to keep.  Leaving blank will import all columns
#' @export
function serial_merge(
	merge_paths,
	sample_key = binfotron::get_default_sample_key(),
	import_clms = NULL
){
	import_df = NULL
	for (import_path in merge_paths){
		this_df = fread(import_path, data.table = F)
		if (!sample_key %in% names(this_df)) stop(paste0("This path does not contain your sample_key, ", sample_key,": ", import_path))
		if (!is.null(import_clms) this_df = this_df[, import_clms[import_clms %in% names(this_df)]]
		if (ncol(this_df) == 1) stop(paste0("This file did not have any data you were looking for: ", import_path))
		if(is.null(import_df)){
			import_df = this_df
		} else {
			import_df = merge(import_df, this_df, by = sample_key, all = T)
		}
	}
	
	if(!is.null(import_clms)){
	  missing_clms = c(import_clms)[!(c(import_clms) %in% names(import_df))]
	  if (length(missing_clms) > 0) stop(paste0(c("Your imported data are missing the following columns:", missing_clms), collapse= "\n"))
	}
}

#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' split_disc_val
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Split a data frame into disc and val df's based on values in a 
#'
#' @param import_df Starting data.frame
#' @param set_clm String to indicate which column contains the sample keys
#' @export
function split_disc_val(
	import_df,
	set_clm = "Set"
){
	if (!all(c("Validation","Discovery") %in% import_df[[set_clm]])) stop(paste0("Your set_clm, ", set_clm, ", does not contain 'Discovery' and 'Validation' strings."))
	val_df = import_df[import_df[[set_clm]] == "Validation", ]
	import_df = import_df[import_df[[set_clm]] == "Discovery", ]
	val_df = val_df[, names(val_df)[names(val_df) != set_clm]]
	import_df = import_df[, names(import_df)[names(import_df) != set_clm]]
	return(list(disc_df=import_df, val_df=val_df))
}