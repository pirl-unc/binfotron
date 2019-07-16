# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# screen_counts
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title screen_counts 
#' 
#' @param my_counts return valid counts
#' 
#' @return Returns couns that are not NA and greater than zero.  Factors are 
#' converted to numerics correctly.
#' 
#' @family diversity
#' 
#' @export
screen_counts = function(my_counts){
  if(class(my_counts) %in% c("character", "factor")){ 
    warning("Your counts are factors not numbers.  This was fixed for diveristy calculations, but this is a huge mistake that should be fixed immediately.  Ask someone if you don't know what this warning is talking about.")
    my_counts = as.numeric(as.character(my_counts))
  } else {
    my_counts = as.numeric(my_counts)
  }
  my_counts = my_counts[!is.na(my_counts)]
  my_counts = my_counts[my_counts>0]
  return(my_counts)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# shannon_entropy
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title shannon_entropy 
#' 
#' @param my_counts vector of postive integers
#' @param should_screen_counts Boolean to indicate if the coutns should be screened
#'   for valid data.  Set to false if the data has already been checked by another function.
#' 
#' @return Returns \code{vegan::diversity(my_counts, index = "shannon")} of the 
#' non-zero counts with NA removed.
#' 
#' @family diversity
#' 
#' @export
shannon_entropy = function(my_counts, should_screen_counts = TRUE)
  #https://stat.ethz.ch/pipermail/r-help/2008-July/167112.html
{
  if (should_screen_counts) my_counts %<>% screen_counts()
  
  if (!checkmate::test_numeric(my_counts, lower = 0, min.len = 2)) return(NA)

  vegan::diversity(my_counts, index = "shannon")
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# inv_simpson
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# https://stat.ethz.ch/pipermail/r-help/2008-July/167112.html
#' @title inv_simpson 
#' 
#' @param my_counts vector of postive integers
#' @param should_screen_counts Boolean to indicate if the coutns should be screened
#'   for valid data.  Set to false if the data has already been checked by another function.
#' 
#' @return Returns \code{vegan::diversity(my_counts, index = "invsimpson")} of the 
#' non-zero counts with NA removed.
#' 
#' @family diversity
#' 
#' @export
inv_simpson = function(my_counts, should_screen_counts = TRUE)
  #https://stat.ethz.ch/pipermail/r-help/2008-July/167112.html
{  
  if (should_screen_counts) my_counts %<>% screen_counts()
  
  if (!checkmate::test_numeric(my_counts, lower = 0, min.len = 2)) return(NA)

  vegan::diversity(my_counts, index = "invsimpson")
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# chao1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title chao1 
#' 
#' @param my_counts vector of postive integers
#' @param should_screen_counts Boolean to indicate if the coutns should be screened
#'   for valid data.  Set to false if the data has already been checked by another function.
#' 
#' @return Returns \code{vegan::estimateR(round(my_counts))["S.chao1"]} of the 
#' non-zero counts with NA removed.
#' 
#' @family diversity
#' 
#' @export
chao1 = function(my_counts, should_screen_counts = TRUE){  
  
  if (should_screen_counts) my_counts %<>% screen_counts()
  
  my_counts = round(my_counts)
  class(my_counts) = "integer"
  if (!checkmate::test_integer(my_counts, lower = 0, min.len = 2)) return(NA)
  
  my_return = suppressWarnings(vegan::estimateR(my_counts)["S.chao1"] %>% as.numeric())
  return(my_return)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# evenness
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title evenness 
#' 
#' @param my_counts vector of postive integers
#' 
#' @return Returns \code{shannon_entropy(my_counts)/log(sum(!is.na(my_counts)))} 
#' of the non-zero counts with NA removed.
#' 
#' @family diversity
#' 
#' @export
evenness = function(my_counts, should_screen_counts = TRUE){  
  
  if (should_screen_counts) my_counts %<>% screen_counts()
  
  if (!checkmate::test_numeric(my_counts, lower = 0, min.len = 2)) return(NA)
  
  
  shannon_entropy(my_counts, should_screen_counts = FALSE)/log(sum(!is.na(my_counts)))
}


# This metric seems backwards to me.  diversity index should mean a higher number means 
# more diversity, but here if one clone account for 50 % of the reads it would be 
# really high.  if 50% are needed it would be really low
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dXX_index
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title dXX_index 
#' 
#' @description
#' From vdjtools:
#' "The estimate equals to 1 - n / N, where n is the minimum number of clonotypes 
#' accounting for at least XX% of the total reads and code N is the total number 
#' of clonotypes in the sample. Computes diversity index that equals to one minus 
#' the minimum fraction of clonotypes accounting for at least 50% of the total 
#' reads."
#' 
#' @param my_counts vector of postive integers
#' @param should_screen_counts Boolean to indicate if the coutns should be screened
#'   for valid data.  Set to false if the data has already been checked by another function.
#'   
#' @return Typically called d50_index. Computes the diversity index dXX, where 
#' XX is a specified fraction.  
#' 
#' @family diversity
#' 
#' @export
dXX_index = function(my_counts, my_fraction, should_screen_counts = TRUE){
  
  if (should_screen_counts) my_counts %<>% screen_counts()
  
  if (!checkmate::test_numeric(my_counts, lower = 0, min.len = 2)) return(NA)
  
  ordered_counts = sort(my_counts[my_counts > 0], decreasing = T)
  min_amount = sum(ordered_counts) * my_fraction
  total_clone_number = length(ordered_counts) # N
  if(total_clone_number > 0){
    min_clone_number = 1 # n
    repeat{
      if(sum(ordered_counts[1:min_clone_number]) >= min_amount){
        break
      }
      min_clone_number = min_clone_number + 1
    }
    return(1-min_clone_number/total_clone_number)
  } else {
    return(NA)
  }
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# has_sufficient_abundance_for_entropy
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title has_sufficient_abundance_for_entropy 
#' 
#' @description
#' Over the range of TCGA TCR entropies (0.975-8.1) we calculated the minimum number 
#' of reads needed to get within the 95th CI of the entropy at 100 billion reads. 
#' See Optimize_Diversity_Metrics_CRSV1371 project file, plot_fraction_entropy.R.
#' 
#' @param my_abundance Integer of the total abundance of the counts
#' @param measured_entropy Measured entropy of the sample.  Value should be close to the 
#' range of (0.975-8.1).
#' 
#' @return Boolean of whether this abundance is enough to just use the measured entropy
#' 
#' @family diversity
#' 
#' @export
#' 
has_sufficient_abundance_for_entropy = function(my_abundance, measured_entropy){
  
  my_exp = 1.37
  min_abundance = (my_exp^((measured_entropy+5)^my_exp)) + 512
  
  if(my_abundance >= min_abundance){
    return(TRUE)
  } else{
    return(FALSE)
  }
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# predict_true_entropy_from_diversity
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title predict_true_entropy_from_diversity 
#' 
#' @description
#' Uses an elastic net model to predict true entropy.  Two models are used to cover
#' different abundnce ranges.  2-1024 is one model 1025-65536 is another.  Over 2^16
#' just returns the shannon entropy, which is pretty accurate over this amount. These
#' models were build around tcga tcr distributions that had an entropy range of 0.975 
#' to 8.1 if corrected entropy is outside of this range then the data should be taken 
#' with caution.
#' 
#' @param my_counts vector of postive integers
#' @param should_screen_counts Boolean to indicate if the coutns should be screened
#'   for valid data.  Set to false if the data has already been checked by another function.
#' @param min_abundance Integer to indicate minimun number of counts needed (ie. \code{sum(my_counts)})
#'   to get a valid prediction.  
#' @param max_abundance Integer over which the function will return shannon entropy.
#' 
#' @return Modeled prediction of diversity 
#' 
#' @family diversity
#' 
#' @export
#' 
predict_true_entropy_from_diversity = function(
  my_abundance,  
  my_richness,
  my_d25,
  my_inv_simpson,
  my_chao1,
  my_shannon_entropy,
  my_evenness,
  min_abundance = 2, 
  should_screen_counts = TRUE
  ){
  
  if (!checkmate::test_numeric(my_abundance, lower = 1, len = 1, any.missing = FALSE)) return(NA) 
  if (!checkmate::test_numeric(my_richness, lower = 2, len = 1, any.missing = FALSE)) return(NA) 
  if (!checkmate::test_numeric(my_d25, lower = 0, upper = 1, len = 1, any.missing = FALSE)) return(NA) 

  if (!checkmate::test_numeric(my_inv_simpson, lower = 1, len = 1, any.missing = FALSE)) return(NA) 
  if (!checkmate::test_numeric(my_chao1, lower = 1, len = 1, any.missing = FALSE)) return(NA) 
  if (!checkmate::test_numeric(my_shannon_entropy, lower = 0, len = 1, any.missing = FALSE)) return(NA) 
  
  if (!checkmate::test_numeric(my_evenness, lower = 0, upper = 1, len = 1, any.missing = FALSE)) return(NA) 
  
  # my_abundance = sum(my_counts)
  # my_shannon_entropy = binfotron::shannon_entropy(my_counts, should_screen_counts = FALSE)
  # 
  # 2^16 is the cutoff on the scale of plot_subsampled_diversity_with_corrections.R 
  if(my_abundance < min_abundance) return(NA)
  
  if(binfotron::has_sufficient_abundance_for_entropy(my_abundance, my_shannon_entropy)) return(my_shannon_entropy)
  
  my_df = data.frame(
    Log2_TRB_Abundance = log2(my_abundance + 1), 
    Log2_TRB_Richness = log2(my_richness + 1), 
    TRB_d25 = my_d25, 
    Log2_TRB_Inv_Simpson = log2(my_inv_simpson + 1),
    Log2_TRB_Chao1 = log2(my_chao1 + 1),
    TRB_Shannon_Entropy = my_shannon_entropy, 
    TRB_Evenness = my_evenness
  )
  
  if(my_abundance <= 1024) {
    model_path = binfotron::get_corrected_entropy_rdata_ab8_1024_ent1_8_path()
  } else {
    model_path = binfotron::get_corrected_entropy_rdata_ab1025_65K_ent1_8_path()
  }
  
  load(model_path, verbose = F)
  feature_names = row.names(optimal_model$beta)
  model_beta = as.matrix(optimal_model$beta)
  final_model_input_features = feature_names[abs(model_beta[[1]]) > 0]
  
  if(any(is.na(my_df[1,final_model_input_features]))) {  # can't model with missing values
    
    return(NA)
    
  } else {
    
    library(glmnet)
    return(predict(optimal_model, as.matrix(my_df[,final_model_input_features]), type = "response") %>% as.numeric())
  }
  
}




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# predict_true_entropy_from_counts
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title predict_true_entropy_from_counts 
#' 
#' @description
#' Uses an elastic net model to predict true entropy.  Two models are used to cover
#' different abundnce ranges.  2-1024 is one model 1025-65536 is another.  Over 2^16
#' just returns the shannon entropy, which is pretty accurate over this amount. These
#' models were build around tcga tcr distributions that had an entropy range of 0.975 
#' to 8.1 if corrected entropy is outside of this range then the data should be taken 
#' with caution.
#' 
#' @param my_counts vector of postive integers
#' @param should_screen_counts Boolean to indicate if the coutns should be screened
#'   for valid data.  Set to false if the data has already been checked by another function.
#' @param min_abundance Integer to indicate minimun number of counts needed (ie. \code{sum(my_counts)})
#'   to get a valid prediction.  
#' @param max_abundance Integer over which the function will return shannon entropy.
#' 
#' @return Modeled prediction of diversity 
#' 
#' @family diversity
#' 
#' @export
#' 
predict_true_entropy_from_counts = function(my_counts,  min_abundance = 2, should_screen_counts = TRUE){
  
  if (should_screen_counts) my_counts %<>% binfotron::screen_counts()
  
  if (!checkmate::test_numeric(my_counts, lower = 0, min.len = 2)) return(NA) 
  
  # add new model
  # return na if abundance is over 4096 and below 2^17
  # return entropy if abundance is over 2^17
  
  my_abundance = sum(my_counts)
  my_shannon_entropy = binfotron::shannon_entropy(my_counts, should_screen_counts = FALSE)
  
  # 2^16 is the cutoff on the scale of plot_subsampled_diversity_with_corrections.R 
  if(my_abundance < min_abundance) return(NA)
  
  if(binfotron::has_sufficient_abundance_for_entropy(my_abundance, my_shannon_entropy)) return(my_shannon_entropy)
  
  my_df = data.frame(
    Log2_TRB_Abundance = log2(my_abundance + 1), 
    Log2_TRB_Richness = log2(length(my_counts) + 1), 
    TRB_d25 = binfotron::dXX_index(my_counts, my_fraction = 0.25, should_screen_counts = FALSE), 
    Log2_TRB_Inv_Simpson = log2(binfotron::inv_simpson(my_counts, should_screen_counts = FALSE) + 1),
    Log2_TRB_Chao1 = log2(binfotron::chao1(my_counts, should_screen_counts = FALSE) + 1),
    TRB_Shannon_Entropy = my_shannon_entropy, 
    TRB_Evenness = binfotron::evenness(my_counts, should_screen_counts = FALSE)
  )
  
  if(my_abundance <= 1024) {
    model_path = binfotron::get_corrected_entropy_rdata_ab8_1024_ent1_8_path()
  } else {
    model_path = binfotron::get_corrected_entropy_rdata_ab1025_65K_ent1_8_path()
  }
  
  load(model_path, verbose = F)
  feature_names = row.names(optimal_model$beta)
  model_beta = as.matrix(optimal_model$beta)
  final_model_input_features = feature_names[abs(model_beta[[1]]) > 0]
  
  if(any(is.na(my_df[1,final_model_input_features]))) {  # can't model with missing values
    
    return(NA)
    
  } else {
    
    library(glmnet)
    return(predict(optimal_model, as.matrix(my_df[,final_model_input_features]), type = "response") %>% as.numeric())
  }
  
}

