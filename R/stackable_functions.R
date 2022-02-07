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



#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' calculate_gene_signatures
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' 
#' @title Calculate gene signatures on my_dt using a gmt file.
#' 
#' @description
#' Calculates gene signatures on my_dt using a gmt file. The method of calculating the 
#' gene signture can be changed using \code{my_fun}
#'   
#' @param gene_cols Character vector indicating the columns that should be used in the 
#'   gene signature. Uses \code{operatable_columns}.  If \code{NULL} if will take all 
#'   numeric columns.
#' @param gene_element For gene names that have a pipe in them, which position should be used ("1|2 etc"). Integer.
#' @param gmt_file_path String path to the gmt file.
#' @param min_genes Integer for the minimum number of genes that have to be found from a
#'   gene signature to keep it.
#' @param col_names Vector of character strings to name the columns that will have this 
#'   operation performed on them.  Uses \code{\link{operatable_columns}}
#' @param my_dt data.table input
#' @param my_fun Function used to combine the genes for each sample. Default is \code{\link{mean}}
#' @param my_summary Character string to change the default comment that will be appended 
#'   to my_dt.
#' @param only_return_signatures Boolean indicating whether or not all of the non-gene columns
#'   should be included in the result.  Gene columns will be dropped.  This is helpful in 
#'   returning categories with the gene signature scores
#' @param summary_output_path If specified this is where the summary data for making this gene 
#'   signature will go.
#' @param sample_key Character string to specify the column that is the sample key. This 
#'   column will not be operated upon. 
#' @param signatures Specifies which gene signtures on the gmt file will be used. If \code{NULL}
#'   all of the signatures wil be used.
#' @param readme_path Optional path to which the comments will be appended.
#' 
#' @section Todos:
#' \itemize{
#'   \item Should convert a data.frame to a data.table and then convert it back for the output.
#' }
#' 
#' @family gene_signature
#' 
#' @export
calculate_gene_signatures = function(
  my_dt = NULL,
  gene_cols = NULL,
  gene_element = 1,
  gmt_file_path,
  min_genes = 1,
  my_fun = median,
  my_summary = "Gene signatures by default are the median of the genes listed in the signature for each sample.",
  only_return_signatures = TRUE,  # as opposed to the other_data_columns (ie non gene columns)
  summary_output_path = NULL,
  sample_key = get_default_sample_key(),
  signatures = NULL,
  readme_path = NULL,
  na_remove = TRUE
){
  assert_data_frame(my_dt)
  
  function_name = "calculating_gene_signatures"
  
  previous_comments = attributes(my_dt)$comments
  
  text_output = make_intro_text(function_name, my_summary)
  
  # init common function procedures.  see pipe source for definitions
  gmt_file = basename(gmt_file_path)
  
  
  gmt_annotation = housekeeping::import_annotation(gmt_file_path)
  text_output = c("Using gmt file: ", gmt_file)
  if(length(gmt_annotation) > 0){
    text_output %<>% c(gmt_annotation, "")
  }
  
  if(grepl(".gmt.txt$", gmt_file_path) | grepl(".gmt$", gmt_file_path)){
    gsc = import_gmt_as_list(gmt_file_path)
  } else if (grepl(".rdata$", gmt_file_path)){
    loaded_gsc_name = load(gmt_file_path, verbose = T)
    gsc = eval(parse(text= loaded_gsc_name))
    gsc = as.list(gsc)
  } else {
    stop("Not sure what type of file this is: gmt file does not end in .gmt.txt, .gmt or .rdata.")
  }
  
  if(signatures %>% is_not_null()){
    gsc = gsc[signatures]
  }
  
  
  
  if(!only_return_signatures){
    other_data_columns = names(my_dt)
    other_data_columns = other_data_columns[other_data_columns %ni% gene_cols]
  }
  
  # gene_dat = fread(input_path) %>% as.data.frame
  if(is.null(gene_cols)){
    gene_cols %<>% operatable_columns(my_dt, acceptable_classes = c("numeric", "integer"), sample_key = sample_key)
  }

  
  start_time = proc.time()[3] #' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  genes_ids = sapply(gene_cols, function(x){strsplit(x, '|', fixed = TRUE)[[1]][gene_element]}) %>% as.character
  gene_indices = match(gene_cols, names(my_dt))
  names(my_dt)[gene_indices] = genes_ids
  
  sig_stat_cols = c("Gene_Signature", "Total_Genes", "Genes_Found", "Percent_Found")
  sig_stats = data.frame(matrix(nrow = length(gsc), ncol = length(sig_stat_cols)))
  names(sig_stats) = sig_stat_cols
  
  immune_sigs = data.frame(matrix(nrow = nrow(my_dt), ncol = (length(gsc) + 1)))
  immune_sigs[,1] = my_dt[,1]
  names(immune_sigs)[1] = sample_key
  drop_col = c()
  gene_sig_names = names(gsc)
  for(my_index in 1:length(gene_sig_names)){
    
    my_gene_set = gene_sig_names[my_index]
    sig_stats$Gene_Signature[my_index] = my_gene_set
    
    my_genes = unlist(gsc[my_gene_set])
    found_ids = my_genes[my_genes %in% genes_ids]
    
    num_total = length(my_genes)
    sig_stats$Total_Genes[my_index] = num_total
    
    num_found = length(found_ids)
    sig_stats$Genes_Found[my_index] = num_found
    
    sig_stats$Percent_Found[my_index] = specify_decimal(100 * num_found/num_total, 1)
    
    #a(paste(my_gene_set, num_total, num_found, percent_found, sep = "\t"))
    if(num_found >= min_genes){ # could put criteria for dropping more gene sigs here
      subdat = my_dt[ , found_ids, with = FALSE, drop = F]
      
      if(na_remove){
        immune_sigs[[my_index + 1]] = apply(subdat, 1, function(x) {my_fun(x[!is.na(x)])})
      } else {
        immune_sigs[[my_index + 1]] = apply(subdat, 1, my_fun)
      }
      names(immune_sigs)[my_index + 1] = my_gene_set
    } else {
      drop_col = c(drop_col, my_gene_set)
    }
  }
  
  # a(paste("For each sample, average the expression values of all of the individual",
  #         "genes in the signature.  That average is the signature expression for",
  #         "that sample.") %>% wrap_sentence)
  if(length(drop_col) > 0){
    text_output %<>% c(paste0("Dropping gene sigs with under ", min_genes," matching gene(s):"))
    text_output %<>% c(paste0(drop_col, collapse = ", "))
    gene_sig_names = gene_sig_names[gene_sig_names %ni% drop_col]
    immune_sigs[, c(names(immune_sigs)[1], gene_sig_names)]
  }
  
  
  if(!is.null(summary_output_path)){
    if(grepl(".csv$", summary_output_path)){
      fwrite(sig_stats, summary_output_path, sep = ",")
    } else {
      fwrite(sig_stats, summary_output_path, sep = "\t")
    }
  }
  
  announce_total_time(function_name, start_time)
  
  immune_sigs %<>% data.table
  
  if(!only_return_signatures){
    immune_sigs = merge(my_dt[ , other_data_columns, with = FALSE, drop = F],  immune_sigs, by = sample_key)
  }
  
  if(is.null(readme_path)){
    attributes(my_dt)$comments = c(previous_comments, text_output, "")
  } else {
    write(text_output, readme_path, append = TRUE)
  } 
  
  return(immune_sigs)
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
#'
#' @return my_dt with zereos in place of NA's
#'  
#' @export
convert_na_to_zeroes = function(
  my_dt,
  my_summary = "Convert NA to zeroes.",
  readme_path = NULL
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
  
  for (i in names(my_dt))
    my_dt[is.na(get(i)), (i):=0]

  if (converted_dt) my_dt %>% as.data.frame()
  
  announce_total_time(function_name, start_time)
  
  if(is.null(readme_path)){
    attributes(my_dt)$comments = c(previous_comments, text_output, "")
  } else {
    write(text_output, readme_path, append = TRUE)
  }
  
  return(my_dt)
}
