#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' regression
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' @title performs the specified regression on the data
#' 
#' @description
#' The purpose of \code{regression} is to perform a regression on the data 
#' across the range of independant and dependant variables provided.  If m
#' 
#' @details
#' This function utilizes one of either \code{\link{glm}} or \code{\link{coxph}} methods.
#' 
#' @param input_dt A data.table that includes all columns of data needed to do the 
#'   analysis:  \code{names(indep_list)}, dep_vars (for glm), \code{names(inclusion_list)},
#'   event_clm & time_clm (for coxph), \code{names(unique(unlist(model_comparison_list))},
#'   and my_grouping.
#' @param indep_list Required named list of column names to use as the independent 
#'   variable. Names of the list will be used to name the output stats.
#'   Example:  \cr 
#'   \code{my_indep_list = list( # no default \cr 
#'     TRA_Chao1 = c("TRA_Chao1"), \cr 
#'     TRB_Chao1 = c("TRB_Chao1"), \cr 
#'     TCR_Chao1 = c("TRA_Chao1", "TRB_Chao1") \cr 
#'   )}
#' @param base_file_name Character string to prefix the names of the output files.
#' @param combined_group_name Character string to call the combinded groups catagory.
#' @param dep_vars Character vector containing the coumn names. Example: \cr 
#'   \code{codemy_dep_vars = c("Age", "SNV_Log2_Neoantigens", "Indel_Log2_Neoantigens")}
#' @param dep_var_families This character vector should contain the names of the 
#'   families to add to the model.  This isn't used for coxph, since the dependent
#'   variable for coxph is survival.  Possible values here should be of the form: 
#'   \code{'Gamma("identity")'} or \code{'gaussian'} and can include anything accepted 
#'   \code{\link{glm}}.
#' @param event_clm For coxph.  The name of the column from which to draw the event 
#'   information.  The column should only contain integers of 1 and 0.  If specified, 
#'   this column needs to be present in \code{input_dt}.
#' @param fdr_by_columns Deprecated. Multiple PValue columns made this overly complicated. Just use \code{\link{binfotron::calc_fdr}} separately.
#' @param fdr_by_columns_for_model_comp Deprecated. Multiple PValue columns made this overly complicated. Just use \code{\link{binfotron::calc_fdr}} separately.
#' @param fdr_method Deprecated. Multiple PValue columns made this overly complicated. Just use \code{\link{binfotron::calc_fdr}} separately.
#' \code{stats::p.adjust.methods}.
#' @param inclusion_list List to specify the samples that should be kept.
#'   For example \code{list(pathology_T_stage = c('T1', 'T2', 'T3'), is_asian = c(TRUE))} 
#'   would drop samples in which the value for the column named 'pathology_T_stage' 
#'   was not either 'T1', 'T2', or 'T3'. Samples must also have 'is_asian' equal 
#'   to \code{TRUE}.  The names of this list should be column names for \code{input_dt}
#' @param model_comparison_list Optional named list of column names that should be used for a 
#'   full and reduced model comparison.  Every group of coulms on this list will be 
#'   run against each dep_vars indep_list combination.  The reduced model will only 
#'   include the items on the list.  The full modle will include the independent 
#'   varible(s) as well.  The names of this list will be what the model comparison 
#'   will be called. The values of this list should be column names in \code{input_dt}.
#'   Example: \cr 
#'   \code{model_comparison_list = list( \cr 
#'     Age = c("Age"), \cr 
#'     Tissue = c("Tissue"), \cr 
#'     Combined = c("Age", "Tissue") \cr 
#'   )}
#' @param model_function A function to return the model.  Important to set 
#'   \code{data = model_dt} in the function.  Do not set glm family.  This will be
#'   added based on \code{dep_var_families}. See examples below.
#'   glm example: \cr 
#'   \code{model_function =function(dep_var = "", indep_vars = "NULL"){ \cr 
#'     paste0("glm(", dep_var, " ~ ",  paste0(indep_vars, collapse = " + "), ", data = model_dt)") \cr   
#'   }} \cr 
#'   coxph example: \cr 
#'   \code{model_function = function(dep_var = "", indep_vars = "NULL"){ \cr 
#'      paste0("coxph(Surv(", time_clm,", ", event_clm, ") ~ ",  paste0(indep_vars, collapse = " + "), ", data = model_dt)") \cr 
#'   }} \cr 
#' @param my_grouping This string is the name of the column you want to use to split 
#'   the data into groups. If specified, this column needs to present in \code{input_dt}.
#' @param output_dir Path to the output directory. The parent directory to the path must exist.
#' @param save_models Boolean on whether you would like to save the models in individual rds files named <base_file_name>_<group_name>_<indep_list_name> .
#' @param sample_clm String to indicate the name of the column for sample names. Only used to output predictions.
#' @param time_clm For coxph.  The name of the column from which to draw the time 
#'   information.  If specified, this column needs to present in \code{input_dt}.
#' @param write_files Boolean on whether you would like to write the output files.
#' 
#' @return List containing several outputs: 
#' \enumerate{
#'   \item stats - data.table with the results of the model output
#'   \item model_comp - data.table with the full vs reduced model comparisons if \code{model_comparison_list} is provided 
#'   \item readme - An output of the comparisons made.
#' }
#'   
#' @section Writes:
#' \itemize{
#'   \item stats file 
#'   \item model_comp file if \code{model_comparison_list} is provided.
#' }
#' 
#' @section Todos:
#' \itemize{
#'   \item Support GAMLSS and ability to determine its own family??
#'   \item Stats output sometimes has blank lines in it
#' }
#' 
#' @section Limitations:
#' \itemize{
#'   \item Haven't fixed to run with ordinals. 
#'   \item Haven't tried models with interactions with it yet
#' }
#' 
#' @family model
#' 
#' @export
regression = function(
  input_dt,
  indep_list,
  base_file_name = "regression_output",
  clear_readme = TRUE,
  combined_group_name = "All",
  dep_vars = NULL,
  dep_var_families = NULL,
  event_clm = "OS_e",
  # fdr_by_columns = NULL,
  # fdr_by_columns_for_model_comp = NULL,
  fdr_method = NULL,
  inclusion_list = list(),
  model_comparison_list = NULL,
  #put warning not to add glm family. we will handle that based on the dependant variable
  model_function = function(dep_var = "", indep_vars = "NULL"){ # for glm
  	paste0("glm(", dep_var, " ~ ", paste0("`", indep_vars, "`", collapse = " + "), ", data = model_dt)")   
  },
  my_grouping = NULL, # "Tissue"
  output_dir = ".",
  sample_clm = get_default_sample_key(),
  save_models = FALSE,
  time_clm = "OS_d",
  write_files = TRUE,
  include_dep_var_in_prediction_name = FALSE
){
  
  library(checkmate)
  
	if (!is.null(fdr_method)){
		if (fdr_method != "none"){
		  stop("No longer supporting FDR correction within this method. Adding multiple PValue columns made this overly complicated. Just use binfotron::calc_fdr separately")
		}
	}
	
  # output paths
  if(!dir.exists(dirname(output_dir))) stop("The parent directory of output_dir does not exist.") # Does dirname(output_dir) exist?
  dir.create(output_dir, showWarnings = F)
  stats_path = file.path(output_dir, paste0(base_file_name, 
                                            "_stats.tsv"))
  if (file.exists(stats_path)) 
    file.remove(stats_path)
  
  model_comparison_path = file.path(output_dir, paste0(base_file_name,"_models.tsv"))
  if (file.exists(model_comparison_path)) 
    file.remove(model_comparison_path)
  
  model_prediction_path = file.path(output_dir, paste0(base_file_name,"_predictions.tsv"))
  if (file.exists(model_prediction_path)) 
    file.remove(model_prediction_path)
  
  readme_path = file.path(output_dir, paste0(base_file_name,"_readme.txt"))
  if (clear_readme) 
    if (file.exists(readme_path)) 
      file.remove(readme_path)
  
  a = function(...) {
    my_output = paste0(...)
    if (!is.null(readme_path)) {
      write(my_output, readme_path, append = TRUE)
    }
    message(my_output)
  }
  my_files = c(stats_path, model_comparison_path, readme_path)
  indep_columns = unique(unlist(indep_list))
  missing_indep = indep_columns[indep_columns %ni% names(input_dt)]
  if (length(missing_indep) > 0) {
    stop(paste0("input_dt is missing indep_vars columns: \n", paste0(missing_indep, sep = "\n", collapse = "")))
  }
  comp_columns = unique(unlist(model_comparison_list))
  missing_comp = comp_columns[comp_columns %ni% names(input_dt)]
  if (length(missing_comp) > 0) {
    stop(paste0("input_dt is missing model_comparison_list columns: \n", paste0(missing_comp, sep = "\n", collapse = "")))
  }
  indep_comp_overlap = intersect(indep_columns, comp_columns)
  if (length(indep_comp_overlap) > 0) {
    stop(paste0("The following columns are in both the model_comparison_list and indep_list: \n", 
                paste0(indep_comp_overlap, sep = "\n", collapse = "")))
  }
  if ("data.table" %ni% class(input_dt)) {
    input_dt %<>% as.data.table()
  }
  check_data_table(input_dt, min.rows = 2, min.cols = 2)
  check_list(indep_list, min.len = 1)
  check_list(model_comparison_list, null.ok = TRUE)
  check_list(inclusion_list, null.ok = TRUE)
  check_character(base_file_name, max.len = 1)
  check_function(model_function)
  column_classes = sapply(1:ncol(input_dt), function(x) {
    class(input_dt[[x]])
  })
  names(column_classes) = names(input_dt)
  output_models = list()
  if (grepl("coxph", model_function())) {
    a("## Running coxph regression")
    library(survival)
    input_dt = input_dt[complete.cases(input_dt[, .SD, .SDcol = c(event_clm, time_clm)]), ]
    is_coxph = TRUE
    dep_vars = "Survival"
    check_character(time_clm, max.len = 1, null.ok = TRUE)
    check_character(event_clm, max.len = 1, null.ok = TRUE)
    if (time_clm %ni% names(input_dt)) {
      stop(paste0("coxph time_clm, ", time_clm, ", can not be found in input_dt."))
    }
    if (class(input_dt[[time_clm]]) %ni% c("integer", "numeric")) {
      stop(paste0("coxph time_clm, ", time_clm, ", on input_dt should be of class integer or numeric."))
    }
    if (any(input_dt[[time_clm]] < 0)) {
      stop(paste0("coxph time_clm, ", time_clm, ", on input_dt should not have values less than 0."))
    }
    if (event_clm %ni% names(input_dt)) {
      stop(paste0("coxph event_clm, ", event_clm, ", can not be found in input_dt."))
    }
    if (any(input_dt[[event_clm]] %ni% 0:1)) {
      stop(paste0("coxph event_clm, ", event_clm, ", can only have values 1 and 0."))
    }
  } else if (grepl("glm", model_function())) {
    a("## Running glm regression")
    is_coxph = FALSE
    missing_dep = dep_vars[dep_vars %ni% names(input_dt)]
    if (length(missing_dep) > 0) {
      stop(paste0("input_dt is missing dep_vars columns: \n", 
                  paste0(missing_dep, sep = "\n", collapse = "")))
    }
    # check for overlap with other columns
    #  when needed the indep_list and dep_vars could be made to be the same, just need to make sure the
    #  dep var removes itself from indep_list items that contain is when it runs.  will do this when it's needed
    indep_dep_overlap = intersect(indep_columns, dep_vars)
    if (length(indep_dep_overlap) > 0) {
      stop(paste0("The following columns are in both the dep_vars and indep_list: \n", 
                  paste0(indep_dep_overlap, sep = "\n", collapse = "")))
    }
    dep_comp_overlap = intersect(dep_vars, comp_columns)
    if (length(dep_comp_overlap) > 0) {
      stop(paste0("The following columns are in both the dep_vars and model_comparison_list: \n", 
                  paste0(dep_comp_overlap, sep = "\n", collapse = "")))
    }
    
    # dep_vars can't be factors
    # dep_vars need to have more than one level
    for (dep_var in dep_vars) {
      if (column_classes[dep_var] %in% c("factor", "character")) 
        stop(paste0("Dependent variable, ", dep_var, ", cannot be a factor or character class."))
      if (length(unique(input_dt[[dep_var]])) < 2) 
        stop(paste0("Dependent variable, ", dep_var, ", only has one value in it's column."))
    }
    
    if (length(dep_var_families) != 1 & (length(dep_var_families) != length(dep_vars))) 
      stop("Length of dep_var_families must be of length 1, which would be applied to all dep_vars OR the same length of dep_vars")
    if (length(dep_vars) > 1) {
      include_dep_var_in_prediction_name = TRUE
    }
    if (length(dep_var_families) == 1) 
      dep_var_families = rep(dep_var_families, length(dep_vars))
  } else {
    stop("Models strings must start with 'glm' or 'coxph'.")
  }
  
  for (indep_column in indep_columns) {
    if (column_classes[indep_column] == "character") {
      warning(paste0("Independent variable, ", indep_column, 
                     ", shouldn't be a character vector. Changing to a factor."))
      input_dt[[indep_column]] = factor(input_dt[[indep_column]])
    }
    if (length(unique(input_dt[[indep_column]])) < 2) 
      warning(paste0("Independent variable, ", indep_column, 
                     ", only has one value in it's column."))
  }
  for (comp_column in comp_columns) {
    if (column_classes[comp_column] == "character") {
      warning(paste0("Model comparison independent variable, ", 
                     comp_column, ", shouldn't be a character vector. Chaning to a factor."))
      input_dt[[comp_column]] = factor(input_dt[[comp_column]])
    }
    if (length(unique(input_dt[[comp_column]])) < 2) 
      warning(paste0("Model comparison independent variable, ", 
                     comp_column, ", only has one value in it's column."))
  }
  
  # check that all indep_list and dep_vars columns  have more than one factor level
  # check mate  dependant var cannot be factors or characters
  # check mate dep var can't only have one level
  input_dt[input_dt == ""] = NA
  if (!is.null(my_grouping)) {
    if (is.na(my_grouping)) {
      my_grouping = NULL
    } else {
      check_character(my_grouping, max.len = 1)
      if (my_grouping %ni% names(input_dt)) {
        stop(paste0("my_grouping, ", my_grouping, ", can not be found in input_dt."))
      }
    }
  }
  if (!is.null(dep_vars) && all(is.na(dep_vars))) 
    dep_vars = NULL

  # need to catch the error first because we still want the warning result
  try_error_warning = function(try_string_eval, my_env) {
    my_value = NA
    a_warn = FALSE
    an_error = FALSE
    warn_msg = NA
    error_msg = NA
    my_result = tryCatch({
      list(return_value = eval(parse(text = try_string_eval), 
                               envir = my_env))
    }, warning = function(w) {
      return(list(return_value = eval(parse(text = try_string_eval), 
                                      envir = my_env), warning = gsub("\n", "", paste(w))))
    }, error = function(e) {
      message(paste0(e))
      my_error = gsub("\n", "", paste(e))
      if (grepl("multi-argument returns are not permitted", 
                my_error)) {
        my_error = "Multiple errors detected."
      }
      return(list(error = my_error))
    })
    if ("return_value" %in% names(my_result)) {
      my_coefs <- coef(my_result$return_value)
      which_are_na <- names(my_coefs)[is.na(my_coefs)]
      if (length(which_are_na) == length(my_coefs)) {
        my_result <- list(error=paste("All coefficients are NA. Reporting as error."))
        my_result$return_value <- NA
      } else if (length(which_are_na) > 0) {
        my_result$warning = paste("NA result from calculation of coefficient(s)", paste(which_are_na, collapse=", "), sep=": ")
      }
      my_value = my_result$return_value
    }
    if ("warning" %in% names(my_result)) {
      a_warn = TRUE
      warn_msg = my_result$warning
    }
    if ("error" %in% names(my_result)) {
      an_error = TRUE
      error_msg = my_result$error
    }
    return(list(
    	return_value = my_value, 
    	warn = a_warn, 
    	warn_msg = warn_msg, 
    	error = an_error, 
    	error_msg = error_msg
    ))
  }
  
  get_glm_stats = function(my_model) {
    coef_mtrx = summary(my_model)$coefficients
    output_dt = rbindlist(lapply(2:nrow(coef_mtrx), function(coef_index) {
    	# coef_index = 1
      output_list = list()
      # summary$my_model can be different than my_model$coefficients as the former gets rid of coef with NA values
      coefficient_name = row.names(coef_mtrx)[coef_index] 
      output_list["Coefficient_Name"] = coefficient_name
      summary_columns = colnames(coef_mtrx)
      # for logtest for glm you'd need to do a full vs reduced model
      pValue_column = summary_columns[grepl("Pr(>|t|)", summary_columns, fixed = T) |
                                        grepl("Pr(>|z|)", summary_columns, fixed = T)]
      output_list["Coef_Wald_PValue"] = coef_mtrx[coef_index, pValue_column]
      output_list["Coef"] = coef_mtrx[coef_index, "Estimate"] %>% as.numeric
      this_coef_data = my_model$data[[coefficient_name]]
      if (this_coef_data %>% is.null) {
        found_it = FALSE
        if (length(my_model$xlevels) > 0) {
          for (x_index in 1:length(my_model$xlevels)) {
            my_var = names(my_model$xlevels)[x_index]
            my_levels = my_model$xlevels[[x_index]]
            for (my_level in my_levels[2:length(my_levels)]) {
              try_me = paste0(my_var, my_level)
              if (try_me == coefficient_name) {
                found_it = TRUE
                this_coef_data = my_model$data[[my_var]]
                output_list["Coefficient_Name"] = paste0(my_var, ": ", my_level, " vs ", my_levels[1])
                break
              }
            }
            if (found_it) 
              break
          }
        } else {
        	# must be a coefficient
          my_var = gsub("TRUE$", "", coefficient_name)
          this_coef_data = my_model$data[[my_var]]
        }
      }
      this_coef_data = this_coef_data[complete.cases(this_coef_data)]
      coef_class = class(this_coef_data)[1]
      if (coef_class == "character") {
        this_coef_data = factor(this_coef_data)
        coef_class = "factor"
      }
      if (coef_class %in% c("numeric", "integer")) {
        my_n = length(this_coef_data) %>% as.character
      } else if (coef_class == "logical") {
        TF_counts = summary(factor(this_coef_data))
        TF_counts = TF_counts[TF_counts %>% order %>% rev]
        # factors 'n' will be expressed as #TRUE:#FALSE
        my_n = paste0(TF_counts, collapse = ":")
      } else if (coef_class == "factor") {
      	# factors 'n' will be expressed as #inCategory:#inControlCategory (ie the first factor level)
        my_n = paste0(summary(this_coef_data)[coef_index], ":", summary(this_coef_data)[[1]])
      } else if (coef_class == "ordered") {
        ordinal_level = substr(coefficient_name, nchar(coefficient_name), nchar(coefficient_name))
        if (ordinal_level == "L") {
          ordinal_level = "1"
        } else if (ordinal_level == "Q") {
          ordinal_level = "2"
        } else if (ordinal_level == "C") {
          ordinal_level = "3"
        }
        ordinal_level = ordinal_level %>% as.numeric
        #first compares level 2 to level 1
        my_n = paste0(summary(this_coef_data)[ordinal_level + 1], ":", summary(this_coef_data)[[ordinal_level]])
      } else {
        stop("coef_class != numeric or integer or logical. need to step through get_stats for this class")
      }
      output_list["N"] = my_n
      if (coef_class %in% c("logical", "numeric", "integer", "factor")) {
        confint_list = tryCatch({
	        suppressMessages(
	            list(df = confint(my_model, level = 0.95))
	          )
	        }, warning = function(w) {
	          suppressWarnings(
	            list(df = confint(my_model, level = 0.95), warning = paste(w))
	          )
	        }, error = function(e) {
	          error_text = gsub("\n", "", paste(e))
	        }
        )
        

        if("df" %in% names(confint_list)){
          output_list["Lower_CI"] = confint_list$df[coefficient_name, "2.5 %"]
          output_list["Upper_CI"] = confint_list$df[coefficient_name, "97.5 %"]
        }
        
        if("warning" %in% names(confint_list)) output_list["CI_Warn"] = paste0(confint_list$warning)
        
      } else {
        warning_msg = paste0(
        	"Confidence intervals have not been tested for the variable type: ", 
        	coef_class, " (", coefficient_name,")"
        )
        output_list["CI_Warn"] = warning_msg
        warning(warning_msg)
      }
      return(output_list)
    }), use.names = T, fill = T)
    return(output_dt)
  }
  get_coxph_stats = function(my_model) {
    if (my_model$nevent > 0) {
      coef_mtrx = summary(my_model)$coefficients
      output_dt = rbindlist(lapply(1:nrow(coef_mtrx), function(coef_index) {
        output_list = list()
        # summary$my_model can be diefferent than my_model$coefficients as the former gets rid of coef with NA values
        coefficient_name = row.names(coef_mtrx)[coef_index]
        output_list["Coefficient_Name"] = coefficient_name
        original_coef_names = names(my_model$assign)
        if (coefficient_name %ni% original_coef_names) {
        	# this is a category and we need to fix the name
          found_it = FALSE
          for (x_index in 1:length(my_model$xlevels)) {
            my_var = names(my_model$xlevels)[x_index]
            my_levels = my_model$xlevels[[x_index]]
            for (my_level in my_levels[2:length(my_levels)]) {
              try_me = paste0(my_var, my_level)
              if (try_me == coefficient_name) {
                found_it = TRUE
                # factors compare to the first factor by default
                output_list["Coefficient_Name"] = paste0(my_var, ": ", my_level, " vs ", my_levels[1])
                break
              }
            }
            if (found_it) 
              break
          }
        }
        summary_columns = colnames(coef_mtrx)
        pValue_column = summary_columns[grepl("Pr(>|t|)", summary_columns, fixed = T) |
                                          grepl("Pr(>|z|)", summary_columns, fixed = T)]
        output_list["Model_LRT_PValue"] = summary(my_model)$logtest['pvalue']
        
        output_list["Coef_Wald_PValue"] = coef_mtrx[coef_index, pValue_column]
        output_list["Coef"] = coef_mtrx[coef_index, "coef"] %>% as.numeric
        
        output_list["Events"] =  paste0(my_model$nevent)# events:non-events
        output_list["Non_Events"] =  paste0(my_model$n - my_model$nevent)# events:non-events
        
        output_list["Hazard_Ratio"] = summary(my_model)$coefficients[coef_index, "exp(coef)"] %>% as.numeric
        
        
        output_list["HR_Lower_CI"] = summary(my_model,conf.int=0.95)$conf.int[coef_index, 3] %>% as.numeric
        output_list["HR_Upper_CI"] = summary(my_model,conf.int=0.95)$conf.int[coef_index, 4] %>% as.numeric
        output_list["HR_SE"] = summary(my_model)$coefficients[coef_index, "se(coef)"] %>% as.numeric
        
        return(output_list)
      }), use.names = T, fill = T)
    } else {
      output_dt = data.table(
      	Model_LRT_PValue = NA, 
      	Coef_Wald_PValue=NA, 
      	Coef=NA, 
      	N=paste0("0:", my_model$n), 
      	Hazard_Ratio=NA, 
      	HR_Lower_CI=NA, 
      	HR_Upper_CI=NA
      )
    }
    return(output_dt)
  }
  if (is_coxph) {
    required_col = unique(c(
    		my_grouping, 
    		names(inclusion_list), 
    		time_clm, 
    		event_clm, 
    		unlist(indep_list), 
    		unlist(model_comparison_list)
    ))
  } else {
  	# first check if all vars are in input_dt
    required_col = unique(c(
    	my_grouping, 
    	names(inclusion_list), 
    	dep_vars, 
    	unlist(indep_list), 
    	unlist(model_comparison_list)
    ))
  }
  if (sum(required_col %ni% names(input_dt)) > 0) {
    missing_col = required_col[required_col %ni% names(input_dt)]
    stop(paste0("The following required columns were not in the input_dt:\n", 
                paste0(paste0("  * ", missing_col), sep = "", collapse = "\n")))
  }
  if (!is.null(inclusion_list) && !is.na(inclusion_list) && 
      (length(inclusion_list) > 0)) {
    a("")
    a("Applying inclusion_list")
    a("Starting input_dt n=", nrow(input_dt))
    for (i_index in 1:length(inclusion_list)) {
      checking_column = names(inclusion_list[i_index])
      if (checking_column %ni% names(input_dt)) 
        stop(paste0("The inclusion list column '", checking_column, 
                    "' was not in your input data.table.  No samples would be included."))
      input_dt = input_dt[input_dt[[checking_column]] %in% 
                            inclusion_list[[i_index]], ]
      a("* After ", checking_column, ": n=", nrow(input_dt))
    }
    a("")
  }
  if (is.null(my_grouping)) {
    my_groups = combined_group_name
  } else {
    my_groups = levels(factor(input_dt[[my_grouping]]))
    if (length(my_groups) != 1) {
      if (!is.na(combined_group_name) | !is.null(combined_group_name)) 
        my_groups = c(combined_group_name, my_groups)
    }
  }
  if (grepl("family", model_function())) {
    warning("You should not include family argument in your model string.  Please put this infomation in dep_var_families.")
  }
  if (!grepl("data = model_dt", model_function(), fixed = T)) {
    stop("Your model string must include 'data = model_dt'")
  }
  get_stats = function(my_model) {
    if (is_coxph) {
      return(get_coxph_stats(my_model))
    } else {
      return(get_glm_stats(my_model))
    }
  }
  pvalue_dt = data.table()
  comp_dt = data.table()
  predictions_dt = data.table()
  show_group_in_model_warning = TRUE
  for (group_index in 1:length(my_groups)) {
  	# group_index = 1
    my_group = my_groups[group_index]
    this_groups_model_comparison_list = list()
    # 0 will be all groups
    group_dt = input_dt
    if (my_group != combined_group_name) {
    	# using get here worked when running the function outside of package but not from inside
    	# if the group is in the model comparison we need to remove it 
      group_dt = group_dt[group_dt[[my_grouping]] == my_group, ]
      if (my_grouping %in% unique(unlist(model_comparison_list))) {
        if (show_group_in_model_warning) {
          warning(paste0("The grouping is included in the model comparisons. ", 
                         "This will crash since it will only have one level when ", 
                         "the individual groups are run. Removing model comparison ", 
                         "that contain the grouping: ", my_grouping))
          show_group_in_model_warning = FALSE # no reason to show this more than once
        }
        doesnt_include_grouping = sapply(
        	model_comparison_list, function(x) { my_grouping %ni% x }
        )
        this_groups_model_comparison_list = model_comparison_list[which(doesnt_include_grouping)]
      }
    } else {
      this_groups_model_comparison_list = model_comparison_list
    }
    a(paste0("- Group: ", my_group, " ------------"))
    for (dep_index in 1:length(dep_vars)) {
    	# dep_index = 1
      dep_var = dep_vars[dep_index]
      a(paste0("  - Dependent variable: ", dep_var))
      if (is_coxph) {
      	# make subdat for this dependent variable
        dep_var_dat = group_dt[complete.cases(group_dt[,.SD, .SDcols = c(time_clm, event_clm)]), ]
      } else {
      	# glm model string
        dep_var_dat = group_dt[complete.cases(group_dt[[dep_var]]), ]
        num_levels = length(levels(factor(as.character(dep_var_dat[[dep_var]]))))
        if (nrow(dep_var_dat) < 2) {
          a(paste0("    - Dependent variable, '", dep_var, "', has fewer than two complete cases. Skipping"))
          next
        } else if (num_levels < 2) {
          a(paste0("    - Dependent variable, '", dep_var, "', has fewer than two levels. Skipping"))
          next
        }
      }
      model_function_w_family = function(raw_string) {
        if (!is_coxph) {
          return(gsub(")$", paste0(", family = ", dep_var_families[dep_index], " )"), raw_string))
        } else {
          return(raw_string)
        }
      }
      for (indep_index in 1:length(indep_list)) {
      	# indep_index = 1
        model_dt = dep_var_dat[
        	complete.cases(
        		dep_var_dat[, unique(c(indep_list[[indep_index]], unlist(this_groups_model_comparison_list))), with = FALSE]
        	), 
        ]
        my_indep_name = names(indep_list)[indep_index]
        a(paste0("    - Independent variable: ", my_indep_name))
        indep_var = indep_list[[indep_index]]
        if (!is_coxph && length(unique(model_dt[[dep_var]])) < 
            2) {
          a(paste0("      - Dependent variable, '", dep_var, "', has fewer than two levels for this independent variable. Skipping"))
          next
        }
        if (nrow(model_dt) < 2) {
          a(paste0("      - Independent variable, '", indep_var, "', has fewer than two complete cases (including model comparisons). Skipping"))
          next
        } else {
          enough_indep_levels = FALSE
          for (this_indep in indep_var) {
            if (length(unique(model_dt[[this_indep]])) > 
                1) {
              enough_indep_levels = TRUE
              break
            }
          }
          if (!enough_indep_levels) {
            a(paste0("      - Independent variable, '", my_indep_name, "', has fewer than two levels for this dependent variable. Skipping"))
            next
          }
        }
        model_dt = droplevels(model_dt)
        model_str = model_function_w_family(model_function(dep_var, 
                                                           indep_var))
        try_model = try_error_warning(model_str, my_env = environment())
        if (try_model$error) {
          my_dt = data.table(
          	Independent = my_indep_name, 
          	Dependent = dep_var, 
          	Group = my_group, 
          	Independent_Var = paste0(indep_var, collapse = ","), 
          	String = model_str, 
          	Error = try_model$error, 
          	Warning = try_model$warn, 
          	Error_Msg = try_model$error_msg, 
          	Warning_Msg = try_model$warn_msg
          )
          a("      * ERROR: ", try_model$error_msg)
        } else {
          my_model = try_model$return_value
          if (include_dep_var_in_prediction_name) {
            model_name = paste0(my_group, "__", dep_var, "_vs_", my_indep_name)
          } else {
            model_name = paste0(my_group, "__", my_indep_name)
          }
          model_name %<>% gsub(" ", "_", .)
          output_models[model_name] = list(my_model)
          if (save_models) 
            saveRDS(model_name, file = file.path(output_dir, paste0(base_file_name, "_", model_name, "_model.rds")))
          
          my_dt = data.table(
          	Independent = my_indep_name, 
          	Dependent = dep_var, 
          	Group = my_group, 
          	Independent_Var = paste0(indep_var, collapse = ","), 
          	get_stats(my_model = my_model), 
          	String = model_str, 
          	Error = try_model$error, 
          	Warning = try_model$warn, 
          	Error_Msg = try_model$error_msg, 
          	Warning_Msg = try_model$warn_msg
          )
          
          prediction_clm = paste0(model_name, "_Prediction")
          model_dt[[prediction_clm]] = predict(my_model, model_dt)
          prediction_model_str = model_function_w_family(model_function(dep_var, prediction_clm))
          try_model = try_error_warning(prediction_model_str, my_env = environment())
          my_model = try_model$return_value
          if (nrow(predictions_dt) == 0) {
            predictions_dt = model_dt[, c(sample_clm, prediction_clm), with = F]
          } else {
          	predictions_dt = merge(
          		predictions_dt, 
          		model_dt[,c(sample_clm, prediction_clm), with = F], 
          		by = sample_clm, all = T
          	)
          }

          # if we predict using the model then we can get the effect size
          if ( length(indep_var) > 1 ){
          	
          	if ( try_model$error ){
          		combined_dt = data.table(
          			Independent = my_indep_name,
          			Dependent = dep_var,
          			Group = my_group,
          			Independent_Var = paste0(indep_var, collapse = ","),
          			Coefficient_Name = prediction_clm,
          			String = prediction_model_str,
          			Error = try_model$error,
          			Warning = try_model$warn,
          			Error_Msg = try_model$error_msg,
          			Warning_Msg = try_model$warn_msg,
          			# complete model because it's on the prediction clm
          			Complete_Model = TRUE 
          		)
          		
          		if ( is_coxph ){
          			combined_dt$Events = sum(model_dt[[event_clm]] == 1)
          			combined_dt$Non_Events = sum(model_dt[[event_clm]] == 0)
          		} else {
          			combined_dt$N = nrow(model_dt)
          		}
          		
          	} else {
          		
          		my_model = try_model$return_value
          	  
          		combined_dt = data.table(
          			Independent = my_indep_name,
          			Dependent = dep_var,
          			Group = my_group,
          			Independent_Var = paste0(indep_var, collapse = ","),
          			get_stats(my_model = my_model),
          			String = prediction_model_str,
          			Error = try_model$error,
          			Warning = try_model$warn,
          			Error_Msg = try_model$error_msg,
          			Warning_Msg = try_model$warn_msg,
          			# complete model because it's on the prediction clm
          			Complete_Model = TRUE 
          		)
          	}
            my_dt$Complete_Model = FALSE # need a way to indicate that the reported pvalue only represent some on the coefficients
            
            my_dt = rbindlist(list(my_dt, combined_dt), use.names = T, fill = T)
            
          } else {
            my_dt$Complete_Model = TRUE # need a way to indicate that the reported pvalue only represent some on the coefficients
          }
          
          
          if(!is.null(this_groups_model_comparison_list) && length(this_groups_model_comparison_list) > 0){
            # for( comp_index in 1:length(this_groups_model_comparison_list)){
            this_comp_dt = rbindlist(lapply(1:length(this_groups_model_comparison_list), function(comp_index){
              my_comp_var = this_groups_model_comparison_list[[comp_index]]
              my_comp_name = names(this_groups_model_comparison_list)[comp_index]
              
              a(paste0("      - Model comparison: ", my_comp_name))
              
              full_model_str = model_function_w_family(model_function(dep_var, c(my_comp_var,indep_var)))
              try_full = try_error_warning(full_model_str, my_env = environment())
              full_model = try_full$return_value
              
              reduced_model_str = model_function_w_family(model_function(dep_var, c(my_comp_var)))
              try_reduced = try_error_warning(reduced_model_str, my_env = environment())
              reduced_model = try_reduced$return_value
              
              if(try_full$error && try_reduced$error){
                return(list(
                  Independent = my_indep_name,
                  Dependent = dep_var,
                  Group = my_group,
                  Model_Comparison = my_comp_name,
                  Independent_Var = paste0(indep_var, collapse = ","),
                  Comparison_Model_Var = paste0(my_comp_var, collapse = ","),
                  N = nrow(model_dt),
                  Full_Model = full_model_str,
                  Reduced_Model = reduced_model_str,
                  Full_Error = try_full$error,
                  Full_Warning = try_full$warn,
                  Reduced_Error = try_reduced$error,
                  Reduced_Warning = try_reduced$warn,
                  Full_Error_Msg = try_full$error_msg,
                  Full_Warning_Msg = try_full$warn_msg,
                  Reduced_Error_Msg = try_reduced$error_msg,
                  Reduced_Warning_Msg = try_reduced$warn_msg
                ))
              } else {
                a1 <- anova(full_model, reduced_model, test='LRT')
                a1_columns = names(a1)
                a1_pValue_column = a1_columns[grepl('Pr(>Chi)', a1_columns, fixed = T) | 
                                                grepl('P(>Chi)', a1_columns, fixed = T) |
                                                grepl('P(>|Chi|)', a1_columns, fixed = T)]
                
                pval <- a1[2, a1_pValue_column]
                logLik_full = logLik(full_model) %>% as.numeric
                logLik_red = logLik(reduced_model) %>% as.numeric # Likelihood is a measure of the extent to which a sample provides support for particular values of a parameter in a parametric model.
                
                return(list(
                  Independent = my_indep_name,
                  Dependent = dep_var,
                  Group = my_group,
                  Model_Comparison = my_comp_name,
                  Independent_Var = paste0(indep_var, collapse = ","),
                  Comparison_Model_Var = paste0(my_comp_var, collapse = ","),
                  Full_LogLik = logLik_full,
                  Reduced_LogLik = logLik_red,
                  Delta_LogLik = logLik_full - logLik_red,
                  Model_LRT_PValue = pval,
                  N = nrow(model_dt),
                  Full_Model = full_model_str,
                  Reduced_Model = reduced_model_str,
                  Full_Error = try_full$error,
                  Full_Warning = try_full$warn,
                  Reduced_Error = try_reduced$error,
                  Reduced_Warning = try_reduced$warn,
                  Full_Error_Msg = try_full$error_msg,
                  Full_Warning_Msg = try_full$warn_msg,
                  Reduced_Error_Msg = try_reduced$error_msg,
                  Reduced_Warning_Msg = try_reduced$warn_msg
                ))

                }
              }), use.names = T, fill = T)
            if (nrow(comp_dt) == 0) {
              comp_dt = this_comp_dt
            } else {
              comp_dt = rbindlist(list(comp_dt, this_comp_dt), 
                                  use.names = T, fill = T)
            }
          }
        }
        if (nrow(pvalue_dt) == 0) {
          pvalue_dt = my_dt
        } else {
          pvalue_dt = rbindlist(list(pvalue_dt, my_dt), 
                                use.names = T, fill = T)
        }
      }
    }
  }

  pvalue_dt$Coefficient_Name = gsub("`", "", pvalue_dt$Coefficient_Name)
  
  pvalue_dt = decode_clms(pvalue_dt, skip_clms = "String")
  
  if(write_files) fwrite(pvalue_dt, stats_path, quote = FALSE, sep = "\t", col.names = TRUE, na = "NA")
  
  if(nrow(comp_dt) > 0){
  	comp_dt = decode_clms(comp_dt)

    a("Done with regression")
    a("")
    readme_content = readLines(readme_path)
    if (write_files) {
      fwrite(comp_dt, model_comparison_path, quote = FALSE, sep = "\t", col.names = TRUE, na = "NA")
      fwrite(predictions_dt, model_prediction_path, quote = FALSE, sep = "\t", col.names = TRUE, na = "NA")
    } else {
      file.remove(readme_path)
    }
    return(list(stats = pvalue_dt, model_comp = comp_dt, 
                readme = readme_content, predictions = predictions_dt))
  } else {
    a("Done with regression")
    a("")
    readme_content = readLines(readme_path)
    if (write_files) {
      fwrite(predictions_dt, model_prediction_path, quote = FALSE, sep = "\t", col.names = TRUE, na = "NA")
    } else {
      file.remove(readme_path)
    }
    return(list(
    	stats = pvalue_dt, 
    	readme = readme_content, 
    	models = output_models, 
    	predictions = predictions_dt
    ))
  }
}
