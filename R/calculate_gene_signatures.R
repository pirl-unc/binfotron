
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
    my_fun = mean,
    my_summary = "Gene signatures by default are the mean of the genes listed in the signature for each sample.",
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
