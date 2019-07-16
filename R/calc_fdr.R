
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' calc_fdr
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' @title Performs an FDR calculation on the fdr_on_columns.
#' 
#' @description
#' The purpose of \code{calc_fdr} is to add FDR correction columns to a data.table using
#' the data.table in data.table out methodology
#' 
#' @param my_dt A data.table or data.frame that includes all columns of data needed to do the 
#'   correction:  \code{fdr_by_columns} & \code{fdr_on_columns}
#' @param fdr_by_columns Optional character vector of column names to specify what 
#' should be group together in the fdr correction using the data.table "by" option.  
#' Leaving this blank will just apply the correction to all of the data values.
#' @param fdr_method String to tell what method to use to FDR correct. Must be one 
#' of the values in \code{stats::p.adjust.methods}.
#' @param fdr_on_columns Character vector to indicate which column(s) should be FDR 
#' corrected.
#' @param readme_path If path is provided the comments from running this function
#' will be appended to the path.
#' 
#' @return data.table or data.frame (depending on the input) with an FDR correction 
#' column for each \code{fdr_on_columns}
#' 
#' @export
calc_fdr = function(
  my_dt,
  fdr_by_columns = NULL,
  fdr_method = "BH",
  fdr_on_columns = "pValue",
  readme_path = NULL
){
  
  # print the text ot the readme and to the console
  a = function(...){
    my_output = paste0(...)
    if(!is.null(readme_path)){
      write(my_output, readme_path, append = TRUE)
    }
    message(my_output)
  }
  
  a("Running calc_fdr")
  a("")
  
  if(fdr_method %in% p.adjust.methods[p.adjust.methods != "none"]){
    a("calc_fdr will correct for False Discovery Rate (FDR) using the stats::p.adjust method, '", fdr_method, "'.")
    a("")
    started_with_data_frame = identical(class(my_dt), "data.frame")
    if(started_with_data_frame){
      my_row_names = row.names(my_dt)
      my_dt %<>% as.data.table()
    }
    
    my_dt$original_row_order = 1:nrow(my_dt)
    
    # make sure all of the fdr_by_columns are there
    if(!is.null(fdr_by_columns)){
      missing_fdr_by_columns = fdr_by_columns[ fdr_by_columns %ni% names(my_dt)]
      if(length(missing_fdr_by_columns)){
        a(paste0(c("The following fdr_by_columns could not be found in input data and will be dropped from the FDR grouping: ", 
                       missing_fdr_by_columns), collapse = "\n  * "))
        
        fdr_by_columns = fdr_by_columns[fdr_by_columns %ni% missing_fdr_by_columns]
      }
      if(length(fdr_by_columns) == 0){
        fdr_by_columns = NULL
      } else {
        a(paste0(c("The FDR corrections will be grouped by: ", 
                   fdr_by_columns), collapse = "\n  * "))
      }
      a("")
    }
    
    for(fdr_on_column in fdr_on_columns){
      fdr_column_is_there = fdr_on_column %in% names(my_dt)
      if(fdr_column_is_there){
        a(paste0("Running on the column, ", fdr_on_column ,"."))
        
        fdr_to_column = paste0(fdr_method, "_", fdr_on_column)
            
        if(is.null(fdr_by_columns)) { 
          # do the fdr correction without grouping the columns
          set(my_dt, j = fdr_to_column, value = p.adjust(p = my_dt[[fdr_on_column]], method = fdr_method))
        } else {
          
          # kill the old column if it's there to make for easier debugging
          if(fdr_to_column %in% names(my_dt)) my_dt[, (fdr_to_column):=NULL]
          
          original_name_order = names(my_dt)
          keep_old_names = names(my_dt)[names(my_dt) %ni% c(fdr_to_column, fdr_by_columns)]
          keep_old_column_text = paste0(names(my_dt), "=", names(my_dt), collapse = ", ")
          fdr_column_text = paste0(fdr_to_column," = p.adjust(",fdr_on_column,", '",fdr_method,"')")
          by_column_text = paste0(fdr_by_columns, collapse = ", ")
          # message(paste0("class(my_dt): ", class(my_dt)))
          parse_text = paste0("my_dt = my_dt[,.( ",
                              keep_old_column_text, ", ", 
                              fdr_column_text,
                              " ), by = .(", by_column_text,")]")
          # message(paste0("parse_text: ", parse_text))
          eval(parse(text = parse_text))
          my_dt = my_dt[, c(original_name_order, fdr_to_column), with = FALSE]
          
        }
        
        # now move the fdr column to right after its pValue column
        #  we know that the new column is last and we know that there is at least 
        #  one column after the pvalue column because we put it there
        pvalue_col_index = which(names(my_dt) == fdr_on_column)
        setcolorder(my_dt, c(1:pvalue_col_index, ncol(my_dt), (pvalue_col_index + 1):(ncol(my_dt)-1)))
      } else {
        a(paste0("calc_fdr will skip applying the FDR correction on the column, ",fdr_on_column ,".  This column was not present in the input data."))
      }  
    }
    my_dt = my_dt[order(my_dt$original_row_order), ]

    # get rid of this extra columns we added to keep the initial row ordering
    my_dt = my_dt[, names(my_dt)[names(my_dt) %ni% "original_row_order"], with = FALSE]

    if(started_with_data_frame){
      my_dt %<>% as.data.frame()
      row.names(my_dt) = my_row_names
    }
  } else {
    a("fdr_method is set to 'none'.  No correction will be performed by calc_fdr.")
  }
  a("")
  a("Done with FDR correction.")
  a("")
  
  return(my_dt)
}

