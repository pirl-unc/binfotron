#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' soft_scale
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' 
#' @title Winsorized scaling
#' 
#' @description
#' Replaces outliers beyond a certain number of standard deviations, ie cap, with a value at that sd.
#'   
#' @param input_data Numeric vector to be scaled
#' @param cap Number indicating the number of sd that you'd like the value to be limited
#' 
#' @export
soft_scale <- function(input_data, cap = scale_cap) {
  # input_data = my_df[[clm]][my_df[[group_clm]] == "Discovery"]
  
  # let's set na's aside from the start
  non_na_indices = which(!is.na(input_data))
  non_na_data = input_data[non_na_indices]
  
  if ( length( unique(non_na_data) ) == 1 ){
    return(rep(0, length(input_data)))
  } else if ( length( unique(non_na_data) ) == 2 ){
    return(as.vector(scale(input_data)))
  } 
  
  # Step 1: Calculate initial z-scores
  initial_zscores <- (non_na_data - mean(non_na_data)) / sd(non_na_data)
  
  # Step 2: Identify outliers (|z| > cap)
  outliers <- abs(initial_zscores) > cap
  
  # Step 3: Calculate mean and SD of non-outliers
  mean_non_outliers = mean(non_na_data[!outliers])
  sd_non_outliers = sd(non_na_data[!outliers])
  
  # Step 4: Replace outliers with capped boundary values
  upper_bound = mean_non_outliers + cap * sd_non_outliers
  lower_bound = mean_non_outliers - cap * sd_non_outliers
  
  adjusted_data = ifelse(non_na_data > upper_bound, upper_bound, ifelse(non_na_data < lower_bound, lower_bound, non_na_data))
  
  # Step 5: Apply scale() on adjusted my_data using the full dataset
  final_scaled_data = input_data
  final_scaled_data[non_na_indices] = as.vector(scale(adjusted_data))
  
  # Return the scaled my_data as a vector
  return(final_scaled_data)
}