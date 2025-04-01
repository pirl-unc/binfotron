# This script defines functions that can be stacked in any order and applied to data.tables
# All functions here should be data.table in and data.table out.  They should also preserve the 
# the comments attribute and add to it to pass along the steps that were done to the data.



#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' get_default_sample_key
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' 
#' @title Return the default f the column name for the default sample key
#' 
#' @description
#' Created to access the default sample key column. unless otherwise specifed this is 
#' the column that will be assumed to contain the sample names
#'
#' @param none
#' 
#' @export
get_default_sample_key = function(){
  return("Run_ID")
  }



#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' add_comments
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' 
#' @title Append my_comment to my_dt 'comments' attribute 
#' 
#' @description
#' Created to attach and build on commnets for data.tables. This allows them to be passed 
#'   along my_dt though the stacks.
#'
#' @param my_comment Vector of character strings 
#' @param my_dt data.table input
#' @param readme_path Optional path to which the comments will be appended.
#' 
#' @export
add_comments = function(
  my_dt = NULL,
  my_comment,
  readme_path = NULL
  ){
  assert_data_table(my_dt)
  
  # cat(paste0(my_comment, "\n"))
  if(is.null(readme_path)){
    attributes(my_dt)$comments = c(attributes(my_dt)$comments, my_comment, "")
  } else {
    write(my_comment, readme_path, append = TRUE)
  }
  
  return(my_dt)
}


#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' scale_columns
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' 
#' @title Scales data.table columns.
#' 
#' @description
#' A stackable wrapper for \code{\link{scale}} to turn it into a stackable transformation
#'
#' @param col_names Vector of character strings to name the columns that will have this 
#'   operation performed on them.  Uses \code{\link{operatable_columns}}
#' @param my_dt data.table input
#' @param my_summary Character string to change the default comment that will be appended 
#'   to my_dt.
#' @param sample_key Character string to specify the column that is the sample key. This 
#'   column will not be scaled.
#' @param should_scale Boolean to indicate whether the data should be scaled or not.
#'   See \code{\link{scale}} 'scale' argument.
#' @param should_center Boolean to indicate whether the data should be centered or not. 
#'   See \code{\link{scale}} 'center' argument.
#' @param readme_path Optional path to which the comments will be appended.
#' 
#' @export
scale_columns = function(
  my_dt = NULL,
  col_names = NULL,
  my_summary = "Scale the data by column.",
  sample_key = get_default_sample_key(),
  should_scale = TRUE,
  should_center = TRUE,
  readme_path = NULL
  ){
  assert_data_frame(my_dt)
  
  function_name = "scale_columns"
  
  previous_comments = attributes(my_dt)$comments
  
  text_output = make_intro_text(function_name, my_summary)
  
  col_names %<>% operatable_columns(my_dt, acceptable_classes = c("numeric", "integer"), sample_key = sample_key)
  
  start_time = proc.time()[3]
  # https://stackoverflow.com/questions/16846380/how-to-apply-same-function-to-every-specified-column-in-a-data-table
  my_dt[ , (col_names) := lapply(.SD, scale, center = should_center, scale = should_scale), .SDcols = col_names]
  #   for (j in col_names) set(my_dt, j = j, value = scale(my_dt[[j]], center = TRUE, scale = FALSE)) # this was much slower: 162.404s
  # More options here: https://stackoverflow.com/questions/16943939/elegantly-assigning-multiple-columns-in-data-table-with-lapply?noredirect=1&lq=1
  # As @Frank pointed to, if you have a lot of columns (say 10,000 or more) then the small overhead 
  # of dispatching to the [.data.table method starts to add up). To eliminate that overhead that there 
  # is data.table::set which under ?set is described as a "loopable" :=. I use a for loop for this type of operation. 
  # It's the fastest way and is fairly easy to write and read.
  # dt[ , names(dt)[20:100] :=lapply(.SD, function(x) sqrt(x) ) , .SDcols=20:100] #  <- good way to use indices 
  # for (col_name in col_names) my_dt[ , (col_name) := scale(my_dt[[col_name]], center = TRUE, scale = FALSE)] # took a loooong time
  announce_total_time(function_name, start_time)
  
  if(is.null(readme_path)){
    attributes(my_dt)$comments = c(previous_comments, text_output, "")
  } else {
    write(text_output, readme_path, append = TRUE)
  }
  
  return(my_dt)
}


#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' round_data
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' 
#' @title Round data.table columns
#' 
#' @description
#' A stackable wrapper for \code{\link{round}} to turn it into a stackable transformation
#'
#' @param col_names Vector of character strings to name the columns that will have this 
#'   operation performed on them.  Uses \code{\link{operatable_columns}}
#' @param digits Integer to specify how mnay decimal digits to round to.
#'   See \code{\link{round}} 'digits' argument.
#' @param my_dt data.table input
#' @param my_summary Character string to change the default comment that will be appended 
#'   to my_dt.
#' @param sample_key Character string to specify the column that is the sample key. This 
#'   column will not be operated upon.
#' @param readme_path Optional path to which the comments will be appended.
#' 
#' @export
round_data = function(
  my_dt = NULL,
  col_names = NULL,
  digits = 0,
  my_summary = "Rounding data.",
  sample_key = get_default_sample_key(),
  readme_path = NULL
  ){ # null will apply to all columns
  assert_data_frame(my_dt)
  
  function_name = "round_data"
  
  previous_comments = attributes(my_dt)$comments
  
  text_output = make_intro_text(function_name, my_summary)
  
  
  col_names %<>% operatable_columns(my_dt, acceptable_classes = c("numeric", "integer"), sample_key = sample_key)

  
  start_time = proc.time()[3]
  # https://stackoverflow.com/questions/16846380/how-to-apply-same-function-to-every-specified-column-in-a-data-table
  my_dt[ , (col_names) := lapply(.SD, round, digits = digits), .SDcols = col_names] # 30sec
  # my_dt[ , (col_names) := round(.SD, digits), .SDcols=col_names] # 151 sec
  
  announce_total_time(function_name, start_time)
  
  if(is.null(readme_path)){
    attributes(my_dt)$comments = c(previous_comments, text_output, "")
  } else {
    write(text_output, readme_path, append = TRUE)
  } 
  
  return(my_dt)
}



#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' require_data_in_samples
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' 
#' @title Drop my_dt columns that aren't present in certain percentage_of_samples according
#'   to requirement function
#' 
#' @description
#' Allows you drop columns that don't meet criteria for X% of samples
#'
#' @param col_names Vector of character strings to name the columns that will have this 
#'   operation performed on them.  Uses \code{\link{operatable_columns}}
#' @param my_dt data.table input
#' @param my_summary Character string to change the default comment that will be appended 
#'   to my_dt.
#' @param percentage_of_samples Percentage of samples that will need to meet the requirement 
#'   to keep the column
#' @param requirement Function by which to judge the column. For example the default 
#'   (\code{function(x){sum(x>0)}}) requires that the value be greater than zero in 75 
#'   percent of the samples.
#' @param sample_key Character string to specify the column that is the sample key. This 
#'   column will not be operated upon.
#' @param readme_path Optional path to which the comments will be appended.
#' 
#' @export
require_data_in_samples = function(
  my_dt = NULL,
  col_names = NULL,
  my_summary = "Require data in smaples",
  percentage_of_samples = 75,
  requirement = function(x){sum(x>0)},
  sample_key = get_default_sample_key(),
  readme_path = NULL
){ # null will apply to all columns
  assert_data_frame(my_dt)
  
  function_name = "require_data_in_samples"
  
  previous_comments = attributes(my_dt)$comments
  
  text_output = make_intro_text(function_name, my_summary)
  
  col_names %<>% operatable_columns(my_dt, sample_key = sample_key)
  
  start_time = proc.time()[3]
  # https://stackoverflow.com/questions/16846380/how-to-apply-same-function-to-every-specified-column-in-a-data-table
  row_count = nrow(my_dt)
  my_test = apply(my_dt[,col_names, with=FALSE], 2, function(x){(100*requirement(x)/row_count) >= percentage_of_samples})
  drop_cols = names(my_test[my_test == FALSE])
  drop_cols = drop_cols[!is.na(drop_cols)]
  if(length(drop_cols) > 0) my_dt[ , (drop_cols) := NULL]
  
  #   for (j in col_names) set(my_dt, j = j, value = scale(my_dt[[j]], center = TRUE, scale = FALSE)) # this was much slower: 162.404s
  # More options here: https://stackoverflow.com/questions/16943939/elegantly-assigning-multiple-columns-in-data-table-with-lapply?noredirect=1&lq=1
  # As @Frank pointed to, if you have a lot of columns (say 10,000 or more) then the small overhead 
  # of dispatching to the [.data.table method starts to add up). To eliminate that overhead that there 
  # is data.table::set which under ?set is described as a "loopable" :=. I use a for loop for this type of operation. 
  # It's the fastest way and is fairly easy to write and read.
  # dt[ , names(dt)[20:100] :=lapply(.SD, function(x) sqrt(x) ) , .SDcols=20:100] #  <- good way to use indices 
  # for (col_name in col_names) my_dt[ , (col_name) := scale(my_dt[[col_name]], center = TRUE, scale = FALSE)] # took a loooong time
  announce_total_time(function_name, start_time)
  
  if(is.null(readme_path)){
    attributes(my_dt)$comments = c(previous_comments, text_output, "")
  } else {
    write(text_output, readme_path, append = TRUE)
  } 
  
  return(my_dt)
}


#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' drop_low_columns
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' 
#' @title Drop my_dt columns are below a certain percentile in value.
#' 
#' @description
#' Specifies a function to calculate the value of the column and a percentile, below which 
#'   samples will be dropped.
#'
#' @param col_names Vector of character strings to name the columns that will have this 
#'   operation performed on them.  Uses \code{\link{operatable_columns}}
#' @param col_value Function that describes how column values are computed
#' @param my_dt data.table input
#' @param my_summary Character string to change the default comment that will be appended 
#'   to my_dt.
#' @param sample_key Character string to specify the column that is the sample key. This 
#'   column will not be operated upon.
#' @param percentile_cutoff Numeric percentile.  Columns below this percentile will be 
#'   dropped
#' @param readme_path Optional path to which the comments will be appended.
#' 
#' @export
drop_low_columns = function(
  my_dt = NULL,
  col_names = NULL,
  col_value = median,
  my_summary = "Drop low expressing columns.",
  sample_key = get_default_sample_key(),
  percentile_cutoff = 10,
  readme_path = NULL
){ # null will apply to all columns
  assert_data_frame(my_dt)
  
  function_name = "drop_low_columns"
  
  previous_comments = attributes(my_dt)$comments
  
  text_output = make_intro_text(function_name, my_summary)
  
  col_names %<>% operatable_columns(my_dt, acceptable_classes = c("numeric", "integer"), sample_key = sample_key)
  
  
  start_time = proc.time()[3]
  my_values = apply(my_dt[,col_names, with=FALSE], 2, function(x){col_value(x)})
  min_value = quantile(my_values, percentile_cutoff/100)
  
  drop_cols = names(my_values[my_values <= min_value])
  my_dt[ , (drop_cols) := NULL]

  announce_total_time(function_name, start_time)
  
  if(is.null(readme_path)){
    attributes(my_dt)$comments = c(previous_comments, text_output, "")
  } else {
    write(text_output, readme_path, append = TRUE)
  } 
  
  return(my_dt)
}


#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' normalize_rows_by_quartile
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' 
#' @title Normalize data values across rows by quartile
#' 
#' @description
#' Normalizes across rows to a certain quartile. Basically wraps \code{\link{quantile}}
#' 
#' @param col_names Vector of character strings to name the columns that will have this 
#'   operation performed on them.  Uses \code{\link{operatable_columns}}
#' @param my_dt data.table input
#' @param my_summary Character string to change the default comment that will be appended 
#'   to my_dt.
#' @param norm_factor Number indicating the value to which the quartile is normalized
#' @param percentile Percentile to which the data will be normalized
#' @param sample_key Character string to specify the column that is the sample key. This 
#'   column will not be operated upon.
#' @param readme_path Optional path to which the comments will be appended.
#' 
#' @export
normalize_rows_by_quartile = function(
  my_dt = NULL,
  col_names = NULL,
  my_summary = "Upper-quartile normalizing data across rows.",
  norm_factor = 1000,
  percentile = 75,
  sample_key = get_default_sample_key(),
  readme_path = NULL
){ # null will apply to all columns
  function_name = "normalize_rows_by_quartile"
  
  
  previous_comments = attributes(my_dt)$comments
  
  text_output = make_intro_text(function_name, my_summary)
  
  col_names %<>% operatable_columns(my_dt, acceptable_classes = c("numeric", "integer"), sample_key = sample_key)
  cat("Found columns\n")

  start_time = proc.time()[3]

  my_dt[,col_names] = t(apply(my_dt[,.SD, .SDcols = col_names], 1, function(x){x * norm_factor/quantile(x[x>0],percentile/100, na.rm=T)})) %>% as.data.table()
  cat(paste0("Time for running ", function_name, ': ', proc.time()[3] - start_time))
  cat("\n\n")
  if(is.null(readme_path)){
    attributes(my_dt)$comments = c(previous_comments, text_output, "")
  } else {
    write(text_output, readme_path, append = TRUE)
  } 
  
  return(my_dt)
}



#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' calc_expression_metrics
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' 
#' @title Quantify expression metrics for qc purposes
#' 
#' @description
#' Calculates normalization factors used for uq normalization.  Also calculates
#' the number of genes that are over \code{min_reads} for QC purposes.
#' 
#' @param col_names Vector of character strings to name the columns that will have this 
#'   operation performed on them.  Uses \code{\link{operatable_columns}}
#' @param my_dt data.table input
#' @param my_summary Character string to change the default comment that will be appended 
#'   to my_dt.
#' @param norm_factor Number indicating the value to which the quartile is normalized
#' @param percentile Percentile to which the data will be normalized
#' @param min_reads Vector of integers for expression to be over to calc percentage of genes over those values
#' @param sample_key Character string to specify the column that is the sample key. This 
#'   column will not be operated upon.
#' @param readme_path Optional path to which the comments will be appended.
#' 
#' @export
calc_expression_metrics = function(
  my_dt = NULL,
  col_names = NULL,
  my_summary = "Calculating expression normalization factors and % of genes with counts over min_reads for QC.",
  norm_factor = 1000,
  percentile = 75,
  min_reads = c(0,1,4),
  sample_key = get_default_sample_key(),
  readme_path = NULL
){ # null will apply to all columns
  function_name = "calc_expression_metrics"
  
  previous_comments = attributes(my_dt)$comments
  
  text_output = make_intro_text(function_name, my_summary)
  
  col_names %<>% operatable_columns(my_dt, acceptable_classes = c("numeric", "integer"), sample_key = sample_key)
  cat("Found columns\n")
  
  start_time = proc.time()[3]
  
  out_dt = my_dt[,sample_key, with = F]
  out_dt$Norm_Factors = apply(my_dt[,.SD, .SDcols = col_names], 1, function(x){norm_factor/quantile(x[x>0],percentile/100, na.rm=T)})
  
  for(min_read in min_reads){
    out_dt[[paste0("Pct_Counts_Over_", min_read)]] = apply(my_dt[,.SD, .SDcols = col_names], 1, function(x){100*sum(x>min_read, na.rm = T)/sum(x>=0, na.rm = T) })
  }
  
  cat(paste0("Time for running ", function_name, ': ', proc.time()[3] - start_time))
  cat("\n\n")
  
  
  if(is.null(readme_path)){
    attributes(out_dt)$comments = c(previous_comments, text_output, "")
  } else {
    write(text_output, readme_path, append = TRUE)
  } 
  
  return(out_dt)
}



#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' log_transform_plus
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' 
#' @title Does a log transformation of the data plus a specified value.
#' 
#' @description
#' This is a wrapper for \code{\link{log}} that also allows the adding of a value prior 
#'   to the transform.
#' 
#' @param add_value Numeric value that will be added prior to the log transform.
#' @param base Numberic value providing the base for the log transform 
#'   See \code{log} 'base.
#' @param col_names Vector of character strings to name the columns that will have this 
#'   operation performed on them.  Uses \code{\link{operatable_columns}}
#' @param my_dt data.table input
#' @param my_summary Character string to change the default comment that will be appended 
#'   to my_dt.
#' @param sample_key Character string to specify the column that is the sample key. This 
#'   column will not be operated upon.
#' @param readme_path Optional path to which the comments will be appended.
#' 
#'   
#' @export
log_transform_plus = function(
  my_dt = NULL,
  add_value = 1,
  base = 2,
  col_names = NULL,
  my_summary = NULL,
  sample_key = get_default_sample_key(),
  readme_path = NULL
){ # null will apply to all columns
  
  function_name = "log_transform_plus"
  
  previous_comments = attributes(my_dt)$comments
  
  text_output = make_intro_text(function_name, my_summary)
  
  
  if(is.null(my_summary)){
    my_summary = paste0("Log", base, " transform (data + ", add_value, ").")
  }

  text_output = make_intro_text(function_name, my_summary)

  col_names %<>% operatable_columns(my_dt, acceptable_classes = c("numeric", "integer"), sample_key = sample_key)
  
  start_time = proc.time()[3]
  
  my_dt[ , (col_names) := lapply(.SD, function(x){log((x + add_value), base = base)}), .SDcols = col_names] # 30sec
  
  announce_total_time(function_name, start_time)
  
  if(is.null(readme_path)){
    attributes(my_dt)$comments = c(previous_comments, text_output, "")
  } else {
    write(text_output, readme_path, append = TRUE)
  } 
  
  return(my_dt) 
}


#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' convert_piped_col_names_to_single_names
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' 
#' @title Convert piped column names into one of the elements of the name.
#' 
#' @description
#' Converts piped column names into one of the elements of the name. 
#'   
#' @param col_names Vector of character strings to name the columns that will have this 
#'   operation performed on them.  Uses \code{\link{operatable_columns}}
#' @param my_dt data.table input
#' @param my_summary Character string to change the default comment that will be appended 
#'   to my_dt.
#' @param position Integer indicating which pipe element do you want
#' @param sample_key Character string to specify the column that is the sample key. This 
#'   column will not be operated upon.
#' @param readme_path Optional path to which the comments will be appended.
#' 
#' @export
convert_piped_col_names_to_single_names = function(
  my_dt = NULL,
  col_names = NULL,
  my_summary = NULL,
  position = 2,
  sample_key = get_default_sample_key(),
  readme_path = NULL
){
  assert_data_frame(my_dt)
  
  function_name = "convert_piped_col_names_to_single_names"
  
  previous_comments = attributes(my_dt)$comments
  
  if(is.null(my_summary)){
    my_summary = paste0("Converting piped column names and returning the name in position ", position, ".")
  }
  
  text_output = make_intro_text(function_name, my_summary)
  
  col_names %<>% operatable_columns(my_dt, sample_key = sample_key)
  
  start_time = proc.time()[3] #' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  col_names = col_names[grepl("|", col_names)] # no point in running this on names that don't have pipes in them
  single_names = sapply(col_names, function(x){strsplit(x, '|', fixed = TRUE)[[1]][position]}) %>% as.character
  col_indices = match(col_names, names(my_dt))
  names(my_dt)[col_indices] = single_names
  
  announce_total_time(function_name, start_time)
  
  if(is.null(readme_path)){
    attributes(my_dt)$comments = c(previous_comments, text_output, "")
  } else {
    write(text_output, readme_path, append = TRUE)
  } 
  
  return(my_dt) 
}



# from https://stackoverflow.com/questions/7235657/fastest-way-to-replace-nas-in-a-large-data-table
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' convert_na_to_zeroes
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' 
#' @title Converts the NA in a data.table or data.frame to zeroes.
#' 
#' @param my_dt Input data.table or data.frame
#' @param my_summary String to describe is doing
#' @param readme_path Path to output readme text to.
#' @param transpose_for_faster_processing Boolean to indicate if you'd like to transpose for potentially faster processing. In a test of a data.frame with 10000 columns and 24 rows where ~75% of values are NA, runtimes were: untransposed: 176 seconds; transposed and transposed back: .05 seconds.  A downside to using this would be if you have different classes in your columns, in which case you could lose formatting or precision.
#'
#' @return my_dt with zereos in place of NA's
#'  
#' @export
convert_na_to_zeroes = function(
		my_dt,
		my_summary = "Convert NA to zeroes.",
		readme_path = NULL,
		transpose_for_faster_processing = F
) {
	
	assert_data_frame(my_dt)
	
	function_name = "convert_na_to_zeroes"
	previous_comments = attributes(my_dt)$comments
	text_output = make_intro_text(function_name, my_summary)
	
	start_time = proc.time()[3] #' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	converted_dt = FALSE
	if(!is.data.table(my_dt)){
		converted_dt = TRUE
		my_dt = as.data.table(my_dt)
	}
	
	if( transpose_for_faster_processing ){
		column_classes <- sapply(my_dt, class)
		
		# Warn if all classes aren't the same.
		if (length(unique(column_classes)) > 1) {
			warning("Not all columns are of the same class. transpose_for_faster_processing may cause issues.")
		}
		
		clm_names = names(my_dt)
		my_dt %<>% t %>% data.table
	}
	
	for (i in names(my_dt))
		my_dt[is.na(get(i)), (i):=0]
	
	if( transpose_for_faster_processing ){
		my_dt %<>% t %>% data.table
		names(my_dt) = clm_names
	}
	
	if (converted_dt) my_dt %>% as.data.frame()
	
	announce_total_time(function_name, start_time)
	
	if(is.null(readme_path)){
		attributes(my_dt)$comments = c(previous_comments, text_output, "")
	} else {
		write(text_output, readme_path, append = TRUE)
	}
	
	return(my_dt)
}


