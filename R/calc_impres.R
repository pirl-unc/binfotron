
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# default_ckpt_pairs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Get the default gene pairs for the impres calculation
#' 
#' @return Returns data frame/table with 2 columns: sample ID & IMPRES score
#' 
#' @export
default_ckpt_pairs = function(){
  gene1_list = c("PDCD1", "CD27", "CTLA4", "CD40", "CD86", "CD28", "CD80", "CD274", "CD86", "CD40",
                 "CD86", "CD40", "CD28", "CD40", "TNFRSF14")
  gene2_list = c("TNFSF4", "PDCD1", "TNFSF4", "CD28", "TNFSF4", "CD86", "TNFSF9", "VSIR", "HAVCR2",
                 "PDCD1", "CD200", "CD80", "CD276", "CD274", "CD86")
  return(as.data.table(cbind(gene1_list, gene2_list), stringsAsFactors = FALSE))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calc_impres
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Calculate an IMPRES score
#' 
#' @param ge_df Gene expression data.frame or data.table. Each column is a gene, each row is a sample.
#'   Must include the sample column indicated by \code{sample_key}.
#' @param id_list List of all sample IDs, list order matches row order in gene expression matrix
#' @param ckpt_pair_df optional data frame with user defined gene comparisons
#' @param require_all_genes Boolean to indicate if the function should proceed if some genes are missing.
#' @param sample_key Character string to specify the column that is the sample key.
#' @param gene_element For gene names that have a pipe in them, which position should be used ("1|2 etc"). Integer.
#' 
#' @return Returns data frame/table with 2 columns: sample ID & IMPRES score
#' 
#' @export
calc_impres = function(
  ge_df,
  ckpt_pair_df = default_ckpt_pairs(),
  gene_element = 1,
  require_all_genes = TRUE,
  sample_key = get_default_sample_key()
){
  
  if (checkmate::checkClass(ge_df, "data.frame") & !checkmate::checkClass(ge_df, "data.table")) {
    input_type = "DF"
    ge_df = as.data.table(ge_df)
    assert(checkmate::checkDataTable(ge_df))
  } else if (checkmate::checkMultiClass(ge_df, c("data.frame", "data.table"))) {
    input_type = "DT"
  } else {
    assert(checkmate::checkDataFrame(ge_df))
  }
  
  if (sample_key %ni% names(ge_df)){
    stop(paste0(sample_key, " must be a column in ge_df."))
  }
  
  if ( "data.frame" %ni% class(ckpt_pair_df) || 
       "data.table" %ni% class(ckpt_pair_df) ){
    stop("ckpt_pair_df must be a data.table or data.frame.")
  }
  ckpt_pair_df = as.data.frame(ckpt_pair_df)

  # get the first part of the names in case these are piped-names 
  names(ge_df) = sapply(names(ge_df), function(x){strsplit(x, '|', fixed = TRUE)[[1]][gene_element]}) %>% as.character()
  
  # check for missing genes first, and list all the ones that are missing
  all_genes = unlist(ckpt_pair_df)
  missing_genes = all_genes[all_genes %ni% names(ge_df)]
  
  if (length(missing_genes) > 0){
    warning(paste0("Your gene data is missing the following columns: \n", paste0(missing_genes, collapse = "\n")))
    if(require_all_genes){
      stop("Terminating calc_impres.")
    } else {
      # drop the missing rows
      cat("Proceeding without the missing gene(s).\n")
      rows_of_missing_genes = unique(c(which(ckpt_pair_df[[1]] == missing_genes), which(ckpt_pair_df[[2]] == missing_genes)))
      ckpt_pair_df = ckpt_pair_df[-rows_of_missing_genes,]
      if (nrow(ckpt_pair_df) == 0) stop("You have no valid gene pair combinations to compute a score.")
    }
  }
  
  
  return_df = ge_df[, sample_key, with = FALSE]
  return_df$IMPRES_Score = 0
  return_df %<>% data.frame()

  for (sample_index in 1:nrow(ge_df)){
    for (pair_index in 1:nrow(ckpt_pair_df)) {
      gene1 = ckpt_pair_df[[1]][pair_index]
      gene2 = ckpt_pair_df[[2]][pair_index]
      gene1_value = ge_df[[gene1]][sample_index]
      gene2_value = ge_df[[gene2]][sample_index]
      if(gene1_value > gene2_value) 
        return_df[sample_index,2] = return_df[[2]][sample_index] + 1
    }
  }
  

  if(input_type == "DF") {
    return(return_df)
  }
  else {
    return(as.data.table(return_df))
  }
}
