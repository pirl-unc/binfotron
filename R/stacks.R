#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' announce_start_of_stack
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Announce the begining of a stacked function
#' @param function_name Character string name of function
#' @export
announce_start_of_stack = function(function_name){return(c("", str_pad(paste0("Start ", function_name, " "), get_annotation_width(), "right", "-")))}

#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' announce_end_of_stack
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Announce the ending of a stacked function
#' @param function_name Character string name of function
#' 
#' @export
announce_end_of_stack = function(function_name){return(c(str_pad(paste0("Finished ", function_name), get_annotation_width(), "left", "-"), ""))}

#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' bgv_lab_prep_count_data_for_gene_signature
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Preps counts data prior to running gene signature.
#' 
#' @description
#' \code{bgv_lab_prep_count_data_for_gene_signature} is the standard Vincent Lab protocol
#'   for preparing gene count data prior to running gene signatures.
#' \itemize{
#'   \item \code{\link{require_data_in_samples}} reuires the gene be expressed in 70% of the samples
#'   \item \code{\link{normalize_rows_by_quartile}} upper quartile normalizes the data 
#'   \item \code{\link{log_transform_plus}} log2 transforms (counts + 1)
#' }
#'   
#' @param gene_cols Character vector indicating the columns that should be used in the gene signature.
#'   Uses \code{operatable_columns}.  If \code{NULL} if will take all numeric columns.
#' @param my_dt data.table input
#'   
#' @family gene_signature
#' 
#' @export
bgv_lab_prep_count_data_for_gene_signature = function(my_dt, gene_cols){
  function_name = "bgv_lab_prep_count_data_for_gene_signature"
  
  my_dt %<>% add_comments(announce_start_of_stack(function_name)) %>%
  
    require_data_in_samples(
      col_names = gene_cols[gene_cols %in% names(.)], 
      my_summary = "Drop gene columns with expression in fewer than 70% of the samples.", # this allows a default setting of the annotation, in case the user doesn't want to
      percentage_of_samples = 70) %>% 
  
    normalize_rows_by_quartile(col_names = gene_cols[gene_cols %in% names(.)]) %>% 
  
    log_transform_plus(col_names = gene_cols[gene_cols %in% names(.)], base = 2, add_value = 1) %>% 

    add_comments(announce_end_of_stack(function_name))

    return(my_dt)
}

#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' bgv_lab_prep_count_data_for_deseq2
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Preps counts data prior to running gene signature.
#' 
#' @description
#' \code{bgv_lab_prep_count_data_for_deseq2} is the standard Vincent Lab protocol
#'   for preparing gene count data prior to running in DESeq2.
#' \itemize{
#'   \item \code{\link{round_data}} Rounds data
#'   \item \code{\link{require_data_in_samples}} Requires genes be in 70% of samples to keep them
#'   \item \code{\link{drop_low_columns log2}} Drops genes that are in the lowest 10% for mean expression.
#' }
#'   
#' @param gene_cols Character vector indicating the columns that should be used in the gene signature.
#'   Uses \code{operatable_columns}.  If \code{NULL} if will take all numeric columns.
#' @param my_dt data.table input
#'   
#' @family differential_expression
#' 
#' @export
#' 
bgv_lab_prep_count_data_for_deseq2 = function(my_dt, gene_cols){
  function_name = "bgv_lab_prep_count_data_for_deseq2"
  
  my_dt %<>% add_comments(announce_start_of_stack(function_name)) %>%
    
    round_data(
      col_names = gene_cols[gene_cols %in% names(.)],
      digits = 0
    ) %>%
    
    # now drop genes that aren't in x number of samples
    require_data_in_samples(
      col_names = gene_cols[gene_cols %in% names(.)], 
      my_summary = "Drop gene columns with expression in fewer than 70% of the samples.", # this allows a default setting of the annotation, in case the user doesn't want to
      sample_key = my_sample_key,
      requirement = function(x){sum(x>0)},
      percentage_of_samples = 70) %>%
    
    # now drop genes with low expression
    drop_low_columns(
      col_names = gene_cols[gene_cols %in% names(.)], my_summary = "Drop columns with medians in the lower 10th percentile.", # this allows a default setting of the annotation, in case the user doesn't want to
      col_value = mean, # this is how column values are computed
      percentile_cutoff = 10) %>% 
    
    add_comments(announce_end_of_stack(function_name))
  
  return(my_dt)
}
