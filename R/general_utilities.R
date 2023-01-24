# stackable_function_dependencies

#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' %ni%
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Not in
#' 
#' @export
`%ni%`<- Negate(`%in%`)


#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' make_intro_text
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Formats intro text for stackable functions
#'  
#' @param function_name Character name of function
#' @param my_summary Summary text to add to function intro text
#'   
#' @export
make_intro_text = function(
  function_name, 
  my_summary = NULL
  ){
  
  text_output = paste0("Running ", function_name, ":")
  text_output = c(text_output, my_summary)
  
  cat("\n")
  cat(paste0(text_output, collapse = "\n"))
  cat("\n")
  
  return(text_output)
}


#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' announce_total_time
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Tells the total time if given a start time.
#' 
#' @param identifier Character name of identifier
#' @param start_time Start time of 
#' 
#' @export
announce_total_time = function(identifier, start_time){
  cat(paste0("Time to run ", identifier, ': ', round(proc.time()[3] - start_time, 2), " seconds"))
  cat("\n\n")
}


#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' operatable_columns
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' 
#' @title Determines what column names can be operated upon.
#' 
#' @description
#' Uses column classes and sample_key to determine what columns can be operated upon.
#'   
#' @param col_names Vector of character strings to name the columns that will have this 
#'   operation performed on them.
#' @param my_dt data.table input
#' @param acceptable_classes Character vector of classes that are acceptable to perform 
#'   whatever operation.  col_names not of this class will not be returned.
#' @param sample_key Character string to specify the column that is the sample key. This 
#'   column will not be operated upon.
#' 
#' @export
operatable_columns = function(
  col_names = NULL,
  my_dt,
  acceptable_classes = NULL,
  sample_key = "Run_ID"
){
  
  if(is.null(col_names)){
    col_names = names(my_dt)
    if(!is.null(sample_key)){
      col_names = col_names[col_names %ni% sample_key]
    }
  }
  
  # ignore columns not present in the data.table
  missing_columns = col_names[ col_names %ni% names( my_dt )]
  num_drop_columns = length(missing_columns)
  
  if(num_drop_columns>0){
    cat(paste0("The following ", num_drop_columns, " columns were not found in the data.table and were droped from subsequent steps:\n"))
    cat(paste0(missing_columns, collapse = ", "))
    cat("\n")
    col_names = col_names[col_names %ni% missing_columns]
  }
  
  # ignore columns that aren't the correct class
  if(!is.null(acceptable_classes)){
    # acceptable_classes = c("integer", "numeric")
    cat(paste0("Requiring column class to be one of the following: ", paste0(acceptable_classes, collapse = ", ")))
    cat("\n")
    
    my_classes = sapply(my_dt[,col_names, with = FALSE], class)
    drop_columns = col_names[my_classes %ni% acceptable_classes]
    num_drop_columns = length(drop_columns)
    if(num_drop_columns > 0){
      cat(paste0("The following ", num_drop_columns, " columns were not the right class and were dropped from subsequent steps:\n"))
      cat(paste0(drop_columns, collapse = ", "))
      cat("\n\n")
      col_names = col_names[my_classes %ni% drop_columns]
    }
  }
  return(col_names)
}




#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' import_gmt_as_list
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Imports a gmt file as a list.
#' 
#' @param gmt_file_path path to gmt file
#' 
#' @export
import_gmt_as_list = function(gmt_file_path){ # outputs as a gene set collection
  
  #get number of lines from the system command
  con <- file(description=gmt_file_path, open="r")
  # com <- paste("wc -l ", gmt_file_path, " | awk '{ print $1 }'", sep='') ]
  # n <- system(command=com, intern=TRUE) %>% as.integer
  
  gmt_import = readLines(con = con, n = -1L, ok = TRUE, warn = FALSE,
                         encoding = "unknown", skipNul = FALSE)
  close(con)
  
  gene_list = list()
  for(line in gmt_import){
    split_line = unlist(strsplit(line,'\t', fixed=T))
    #gene_set_name = format_gene_set_name(split_line[1]) # get rid of periods and hyphens
    gene_set_name = format_gene_set_name(split_line[1]) # get rid of periods and hyphens
    
    gene_set_genes = split_line[3:length(split_line)]
    gene_set_genes = gene_set_genes[gene_set_genes != 'NA']
    gene_set_genes = gene_set_genes[gene_set_genes != '']
    assembled_code = paste0("gene_list$`", gene_set_name, "` = c('", paste(gene_set_genes, collapse="','"), "')")
    #print(line)
    eval(parse(text = assembled_code))
  }
  
  # testList <- list('hsa-mir-451'=c('SATB2', 'MECP2', 'CTNNBIP1', 'SATB2'), 'hsa-mir-452'=c('SATB2', 'MEIS2', 'PRDM16', 'PRDM16'), 'hsa-mir-453'=c('SATB2', 'SNAI1', 'MECP2'))
  
  n <- names(gene_list)
  uniqueList <- lapply(gene_list, unique)# need unique values in list elements to make genesets
  
  #   #make a function to create the sets
  #   makeSet <- function(geneIds, n) {
  #     GeneSet(geneIds, geneIdType=SymbolIdentifier(), setName=n)
  #   }
  #   
  #   #apply the function to each element in the list and make a list of genesets
  #   gsList <- gsc <- mapply(makeSet, uniqueList[], n)
  #   
  #   #make the geneset collection
  #   gsc <- GeneSetCollection(gsList)
  
  return(uniqueList)
  
}



#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' format_gene_set_name
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Swaps hyphens for underscores
#' 
#' @param my_string String from which to swap hyphens
#' 
#' @export
format_gene_set_name = function(my_string) {
  #x <- gsub('\\.|\\-','_', my_string) # swap hyphens for underscores
  x <- gsub('\\-','_', my_string) # swap hyphens for underscores
  return(x)
}



# need to grab all the ones in the folder 
# update_gene_signatures = function() { sub(gene_sig_file_tag, '', update_gene_sig_files()) }
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get_ipres_gene_sigs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Returns the names of the IPRES gene signatures.
#' 
#' @export
get_ipres_gene_sigs = function(){
  return(c("Hugo_Responder","Hugo_FDR_Responder","Vincent_IPRES_Responder",
           
           "Hugo_NonResponder","Hugo_FDR_NonResponder","Vincent_IPRES_NonResponder",
           "Hugo_IPRES26","Hugo_IPRES22","Hugo_IPRES08","Hugo_IPRES06"))
}




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get_gene_sigs_clustered_by_category
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Gets the name of the immune gene signatures ordered by cell types.
#' 
#' @export
get_gene_sigs_clustered_by_category = function(){
  return(
    
    c(
      "KardosChai_ImSuppress","IglesiaVincent_BCell","Palmer_BCell","Bindea_BCells",
      "Schmidt_BCell","GO_BCR_Signaling","Fan_IGG","Rody_TNBC_BCell",
      
      "Palmer_CD8","IglesiaVincent_CD8","Bindea_CD8_TCells","IglesiaVincent_TCell",
      "Palmer_TCell","Bindea_TCells","GO_TCR_Signaling","Rody_TNBC_TCell","Bindea_THelper",
      "Bindea_Th1_Cells","Bindea_Th2_Cells","Bindea_Th17_Cells","Bindea_TReg","Bindea_Cytotoxic_Cells",
      
      "Rody_IL8","IglesiaVincent_CD68","Rody_LCK","Beck_Mac_CSF1","Bindea_Macrophages",
      "IglesiaVincent_MacTh1",
      
      "Prat_Claudin","KardosChai_EMT_DOWN","KardosChai_EMT_UP", "Chan_TIC"
      
      ,"Bindea_aDC","Bindea_DC","Bindea_Eosinophils","Bindea_iDC",
      "Bindea_Mast_Cells","Bindea_Neutrophils",
      "Bindea_NK_CD56bright","Bindea_NK_CD56dim","Bindea_NK_Cells",
      "Bindea_pDC","Bindea_Tcm","Bindea_Tem","Bindea_TFH","Bindea_Tgd",
      
      get_ipres_gene_sigs()
    )
  )
}




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get_variable_groupings
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Lists the immune gene signatures grouped by cell type.
#' 
#' @export
immune_gene_set_groupings = function(){
  return(list(Immunosuppresssion = 'KardosChai_ImSuppress',
              B_Cells     = c("IglesiaVincent_BCell","Palmer_BCell","Bindea_BCells",
                              "Schmidt_BCell","GO_BCR_Signaling","Fan_IGG","Rody_TNBC_BCell"),
              T_Cells     = c("Palmer_CD8","IglesiaVincent_CD8","Bindea_CD8_TCells","IglesiaVincent_TCell",
                              "Palmer_TCell","Bindea_TCells","GO_TCR_Signaling","Rody_TNBC_TCell","Bindea_THelper",
                              "Bindea_Th1_Cells","Bindea_Th2_Cells","Bindea_Th17_Cells","Bindea_TReg","Bindea_Cytotoxic_Cells"),
              Macrophages = c( "Rody_IL8","IglesiaVincent_CD68","Rody_LCK","Beck_Mac_CSF1","Bindea_Macrophages"),
              
              Claudin     = 'Prat_Claudin', 
              EMT_Down    = 'KardosChai_EMT_DOWN', 
              EMT_Up      = 'KardosChai_EMT_UP', 
              Dendritic_cells = c("Bindea_aDC","Bindea_DC","Bindea_Eosinophils","Bindea_iDC"),
              Eosinophils =  'Eosinophils',
              Mast_Cells  = 'Mast_Cells', 
              Neutrophils = 'Neutrophils', 
              NK_Cells    = c('NK_CD56bright_Cells', 'NK_CD56dim_Cells','NK_Cells'),
              IPRES_Responders = c("Hugo_Responder","Hugo_FDR_Responder","Vincent_IPRES_Responder"), 
              IPRES_NonResponders = c("Hugo_NonResponder","Hugo_FDR_NonResponder","Vincent_IPRES_NonResponder",
                                      "Hugo_IPRES26","Hugo_IPRES22","Hugo_IPRES08","Hugo_IPRES06"),
              Other  = c('Eosinophils', "Mast_Cells","Prat_Claudin", "KardosChai_EMT_DOWN","KardosChai_EMT_UP",'Tcm', 'Tem', 'TFH', 'Tgd', 'TIC')
  ))
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pvalue_stars - took from gtools::stars.pval since this isnt' available for 3.3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Turns pvalue into stars 
#' 
#' @description 
#' Modified from \code{\link{gtools::stars.pval}}. 
#' 
#' @param p.value decimal pvalue
#' @param barely_sig decimal pvalues under this return *
#' @param decent_sig decimal; pValues under this return **
#' @param whoa_nelly decimal; pValues below this return ***
#' @param meh decimal; pValues below this return ***
#' 
#' @return pvalue stars
#' \itemize{
#'   \item *** for 0-000 to 0.001
#'   \item **  for 0.001 to 0.01
#'   \item *   for 0.01  to 0.05
#'   \item "" for over 0.05
#' }
#' 
#' @export
pvalue_stars = function(
  p.value, 
  barely_sig = 0.05, 
  decent_sig = 0.01, 
  whoa_nelly = 0.001, 
  bluh = ""){
  if (p.value <= whoa_nelly){
    return("***")
  } else if (p.value <= decent_sig){
    return("**")
  } else if (p.value <= barely_sig){
    return("*")
  } else {
    return(bluh)
  }
  # unclass(symnum(p.value, corr = FALSE, na = FALSE, 
  #                cutpoints = c(0, 0.001, 0.01, 0.05), 
  #                symbols = c("***", "**",  "*")))
}



#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' is_not_null
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Is not null
#' 
#' @description 
#' Negates is.null.  Just reads better than \code{\link{!is.null}}.  
#' 
#' @param my_vector A vector
#' 
#' @return TRUE or FALSE for each item in the vector
#' 
#' @export
is_not_null = function(my_vector) {
  !is.null(my_vector)
}


# https://stackoverflow.com/questions/3443687/formatting-decimal-places-in-r
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' specify_decimal
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title String rounded to a specific decimal places
#' 
#' @param my_number Number that will be rounded
#' @param number_of_digits Integer of didgits to round to.
#' 
#' @export
specify_decimal = function(my_number, number_of_digits=0){
	my_number = as.numeric(my_number)
	format(round(my_number, number_of_digits), nsmall = number_of_digits)
}


# https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' gg_color_hue
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Get ggplot colors
#' 
#' @param n Integer for the number of colors to return
#' 
#' @export
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write_metadata
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Adds metadata to pdf
#' 
#' @param save_path Character string for directory
#' @param pdf_metadata Character vector of pdf metadata to add to pdf
#' 
#' @export
write_metadata = function(
  save_path, 
  pdf_metadata
){
  if (system("which pdftk", intern = T) != 1) {
    cat("Adding metadata to pdf.\n")
    temp_file = strsplit(save_path, "/") %>% unlist
    temp_file[length(temp_file)] = paste0("TEMP_", temp_file[length(temp_file)])
    temp_file = paste0(temp_file, collapse = "/")
    old_meta = system(paste0("pdftk '", save_path, "' dump_data"), 
                      intern = T)
    new_meta = NULL
    for (this_arg in pdf_metadata) {
      this_arg = strsplit(this_arg, ": ") %>% unlist
      this_arg_meta = c("InfoBegin", paste0("InfoKey: ", 
                                            this_arg[1]), paste0("InfoValue: ", this_arg[2]))
      new_meta = c(new_meta, this_arg_meta)
    }
    new_meta = c(new_meta, old_meta)
    temp_file = strsplit(save_path, "/") %>% unlist
    temp_file[length(temp_file)] = paste0("TEMP_", temp_file[length(temp_file)])
    temp_file = paste0(temp_file, collapse = "/")
    temp_txt = strsplit(temp_file, ".", fixed = T) %>% unlist
    temp_txt = paste0(temp_txt[1], ".txt")
    write.table(new_meta, temp_txt, quote = FALSE, sep = "\t", 
                col.names = FALSE, row.names = FALSE)
    system(paste0("pdftk '", save_path, "' update_info '", temp_txt, 
                  "' output '", temp_file, "' dont_ask"))
    system(paste0("mv '", temp_file, "' '", save_path, "'"))
    system(paste0("rm '", temp_txt,"'"))
  }
  else {
    print("Cannot add metadata. pdftk is not installed on this image.")
  }
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# minimal_decimal_places_ordered
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Reduce decimal places needed
#' 
#' @description
#' Reduce the number of decimals down to the minimal needed to preserve the same number of levels in the original data. Returns data as factors to represent values in the same order they were given. 
#' This is handy for whne you 
#'  
#' @param numeric_vector Values to represent
#' @param min_decimal Integer to indicate the number of folds in which to split the data
#' 
#' @export
minimal_decimal_places_ordered = function(numeric_vector, min_decimal = 0){
	# min_decimal = 0
	initial_levels =  unique(numeric_vector)
	initial_level_count = length(initial_levels)
	test_levels = unique(specify_decimal(initial_levels, min_decimal))
	test_level_count = length(test_levels)

	while(initial_level_count > test_level_count){
		min_decimal = min_decimal + 1
		test_levels = unique(specify_decimal(initial_levels, min_decimal))
		test_level_count = length(test_levels)
	}
	numeric_vector %<>% factor(levels = sort(unique(.)))
	levels(numeric_vector) %<>% specify_decimal(min_decimal)
	return(numeric_vector)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# balanced_folds
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Get indices of a fold of data that balances the distribution of values.
#' 
#' @description
#' This was written by Joel Parker.
#'  
#' @param y Vector to balance
#' @param nfolds Integer to indicate the number of folds in which to split the data
#' 
#' @export
balanced_folds<- function(y, nfolds = 3) { # from Joel
	permuteRows = function(x) {
		dd = dim(x)
		n = dd[1]
		p = dd[2]
		mm = runif(length(x)) + rep(seq(n) * 10, rep(p, n))
		matrix(t(x)[order(mm)], n, p, byrow = T)
	}

	totals = table(y)
	fmax = max(totals)
	nfolds = min(nfolds, fmax)
	folds = as.list(seq(nfolds))
	yids = split(seq(y), y)
	bigmat = matrix(NA, ceiling(fmax / nfolds) * nfolds, length(totals))
	for(i in seq(totals))
		bigmat[seq(totals[i]), i] = sample(yids[[i]])
	smallmat = matrix(bigmat, nrow = nfolds)
	smallmat = permuteRows(t(smallmat))
	res = vector("list", nfolds)
	for(j in 1:nfolds) {
		jj = !is.na(smallmat[ , j])
		res[[j]] = smallmat[jj, j]
	}

	return(res)
}

# how is the t.test method different than the z test
# https://sphweb.bumc.bu.edu/otlt/MPH-Modules/PH717-QuantCore/PH717-Module6-RandomError/PH717-Module6-RandomError11.html
# https://www.statskingdom.com/doc_t_test_z_test.html
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# t_test_ci -  calculates mean & CI for residuals using t-distribution
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title calculates mean & CI for residuals using t-distribution
#' 
#' @description
#' This CI method does not assume we know the sd and uses different t-distributions to get a better estimate of the distribution of the data.
#' Both t_test_ci and z_test_ci assume the data are normally distributed. It will usually be better to use the t_test_ci. This will equal the z_test_ci as N -> inf
#' 
#' @param x Numeric vector to correlate with vector_b
#' 
#' @export
t_test_ci = function(x){
  my_mean = mean(x)
	my_result = tryCatch({
		my_ci = t.test(x)$conf.int
		my_result = list(mean = my_mean, lower_ci = my_ci[1], upper_ci = my_ci[2])
	}, error = function(e) {
		my_mean = mean(x)
		my_result = list(mean = my_mean, lower_ci = NA, upper_ci = NA)
	})
	return(my_result)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# z_test_ci -  calculates mean & CI for residuals using z-score sd
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title calculates mean & CI for residuals using z-score sd
#' 
#' @description
#' The z test assumes we know the sd ( ie we have a large number of samples n > 30). Both t_test_ci and z_test_ci assume the data are normally distributed. It will usually be better to use the t_test_ci. This will equal the t_test_ci as N -> inf
#'  
#' @param x Numeric vector to correlate with vector_b
#' @param show_warning Boolean to indicate whether a warning be shown if there are too few samples to accurately use this method to get confidence intervals
#' 
#' @export
z_test_ci = function(x, show_warning = T) {
  if (show_warning) {
    if (length(x) < 30) {
      warning("N is fewer than 30 values so a t-test would be more appropriate in this case. Try using t_test_ci")
    }
  }
  mean = mean(x)
  my_sd = sd(x)/sqrt(length(x))
  lower_ci = (mean - (1.96 * my_sd))
  upper_ci = (mean + (1.96 * my_sd))
  return(list(mean = mean, lower_ci = lower_ci, upper_ci = upper_ci))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get_rho -  calculates Rho & CI
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Run cor.test and passes results 
#'  
#' @param vector_a Numeric vector to correlate with vector_b
#' @param vector_b Numeric vector to correlate with vector_a
#' @param my_method string to pass to cor.test: 'pearson' or 'spearman'
#' 
#' @export
get_rho = function(vector_a, vector_b, my_method = "pearson") {
  test_result = cor.test(vector_a, vector_b,
                         alternative = "two.sided", method = my_method, exact = F)
  conf_int = test_result$conf.int
  my_rho = as.numeric(test_result$estimate)
  return(list(rho = my_rho, lower_ci = conf_int[1], upper_ci = conf_int[2], pvalue = test_result$p.value))
}
  
