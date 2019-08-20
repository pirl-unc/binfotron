# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calc_IMPRES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title calc_IMPRES 
#' 
#' @param func_GE_DF Gene expression matrix, each column is a gene, each row is a sample
#' @param ID_list List of all sample IDs, list order matches row order in gene expression matrix
#' @param user_ckpt_pair_df optional data frame with user defined gene comparisons
#' 
#' @return Returns data frame/table with 2 columns: sample ID & IMPRES score
#' 
#' @export
calc_IMPRES = function(func_GE_df, ID_list, user_ckpt_pair_df) {
  if(checkmate::checkClass(func_GE_df, "data.frame") == "TRUE" & checkmate::checkClass(func_GE_df, "data.table") != "TRUE") {
    input_type = "DF"
    func_GE_df = as.data.table(func_GE_df)
    checkmate::assert(checkmate::checkDataTable(func_GE_df))
  }
  else if(checkmate::checkMultiClass(func_GE_df, c("data.frame", "data.table")) == "TRUE") {
    input_type = "DT"
  }
  else {
    assert(checkmate::checkDataFrame(func_GE_df))
  }
   
  assert(checkmate::checkList(ID_list))
  
  if(!missing(user_ckpt_pair_df)) {
    if(checkmate::checkMultiClass(user_ckpt_pair_df, "data.table") == "FALSE") {
      if (checkmate::checkClass(user_ckpt_pair_df, "data.frame") == "TRUE") {
        ckpt_pair_df = as.data.table(user_ckpt_pair_df)
        colnames(ckpt_pair_df) = c("gene1_list", "gene2_list")
      }
    }
    else {
      checkmate::assert(checkmate::checkDataFrame(user_ckpt_pair_df))
    }
  }
  else {
    gene1_list = c("PDCD1", "CD27", "CTLA4", "CD40", "CD86", "CD28", "CD80", "CD274", "CD86", "CD40",
                   "CD86", "CD40", "CD28", "CD40", "TNFRSF14")
    gene2_list = c("TNFSF4", "PDCD1", "TNFSF4", "CD28", "TNFSF4", "CD86", "TNFSF9", "C10orf54", "HAVCR2",
                   "PDCD1", "CD200", "CD80", "CD276", "CD274", "CD86")
    ckpt_pair_df = as.data.table(cbind(gene1_list, gene2_list), stringsAsFactors = FALSE)
  }
  
  samp_num = nrow(func_GE_df)
  IMPRES_res = vector()

  for (k in 1:samp_num){
    IMPRES = 0
    for (i in 1:nrow(ckpt_pair_df)) {
      tryCatch(
        if (func_GE_df[k, ckpt_pair_df[i,]$gene1, with = FALSE] > func_GE_df[k, ckpt_pair_df[i,]$gene2, with = FALSE]) IMPRES = IMPRES + 1 else IMPRES = IMPRES,
        error=function(warning_message) {
          message("Either ", paste0(ckpt_pair_df[i,]$gene1, " or ", ckpt_pair_df[i,]$gene2, " does not exist in your gene expression matrix. Please check for possible aliases & try again."))
          stop(call. = FALSE)
        }
      )
    }
    IMPRES_res[k] = IMPRES
  }
  
  IMPRES_res_df = as.data.frame(cbind(as.character(ID_list), as.numeric(IMPRES_res)))
  colnames(IMPRES_res_df) = c("Patient", "IMPRES")

  if(input_type == "DF") {
    return(IMPRES_res_df)
  }
  else {
    return(as.data.table(IMPRES_res_df))
  }
}
