#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' model_differential_expression
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' @title Does a differential gene expression analysis. 
#' 
#' @description
#' The purpose of \code{model_differential_expression} is to compare the raw read counts of gene expression data 
#' between different groups of samples to see if there is differential gene expression. This is different from model_gene_expression in that this one 
#' tries to incorporate the changes needed for the LCCC-Bioinformatics group to use. Dumping sam and edger here. Letting
#' old Deseq (ie not deseq2) go for a bit until it's needed.
#' 
#' @details
#' This function utilizes one of either \code{\link{DESeq}} or \code{\link{DESeq2}} methods.  \code{\link{DESeq}} is recommended for single cell data.
#' 
#' @param analysis_method String option out of the following choices:  DESeq2, DESeq, SAM, or edgeR. Indicating which method
#'   should be used to do the analysis.
#' @param base_file_name String to specify the file name.
#' @param core_number Integer to indicate the numboer of cores that should be used.
#' @param deseq2_results_cooksCutoff Set to \code{Inf} or \code{FALSE} to disable the resetting of p-values to \code{NA}. 
#'   Gets passed to \code{\link{DESeq2::results}} 'independentFiltering' argument
#' @param deseq2_results_independentFiltering Gets passed to \code{\link{DESeq2::results}} 'cooksCutoff' argument
#' @param design_var This needs to be setup to allow more complex comparisons and and matching for samples.
#' @param gene_expression_cols Character vector with the names of the columns with genes in them.
#' @param gmt_file_fdr_cutoffs Numeric vector of cutoffs to use for the FDR significant values.  Two gene signtures
#'   will be made of all the genes that have under the fdr pValue: one for up genes and one for down.
#' @param gmt_file_pvalue_cutoffs Numeric vector of cutoffs to use for the pValue significant values.  Two gene signtures
#'   will be made of all the genes that have under the pValue: one for up genes and one for down.
#' @param gmt_ref String indicating what should go in the reference part of the gmt file.
#' @param imported_annotation Character vector to include what steps were done to the data prior to this analysis.  This module will 
#'   add on to those steps.
#' @param my_grouping This string is the name of the column you want to use to split the data into groups.
#' @param output_dir Path to the output directory.  This will be calculated automatically if left blank.  
#' @param sample_key_col String matching the name of the column that should be used to identify the unique sample identifiers.
#' 
#' @return List containing several outputs contains: 
#' \enumerate{
#'   \item heatmap_pack - to be used by \code{\link{view_heatmap}}; 
#'   \item volcano_pack - to be used by \code{\link{view_volcano_plot}}; 
#'   \item gmt_id_path - path to a gmt file with genes represetned by ids (typically entrez);
#'   \item gmt_symbol_path - path to a gmt file with genes represented by symbols (typically hgnc);
#'   \item stats_path - path to stats of Gene_Name, Fold_Change, pValue, FDR_pValue
#' }
#'   
#' @section Writes:
#' \itemize{
#'   \item stats file 
#' }
#' 
#' @section Todos:
#' \itemize{
#'   \item Fix to match by patient ids.
#'   \item Setup design_var and contrasts
#' }
#' 
#' @section Limitations:
#' \itemize{
#'   \item Can only match two samples at this point.
#' }
#' 
#' @section See also:
#' \itemize{
#'   \item \code{\link{view_heatmap}}
#'   \item \code{\link{view_volcano_plot}}
#' }
#' 
#' @family differential_expression, model
#' 
#' @export
model_differential_expression = function(
  my_dt = NULL,
  analysis_method = "DESeq2",
  base_file_name = NULL,
  base_title = NULL,
  core_number = round(detectCores()/2),
  deseq2_results_cooksCutoff = NULL, # Inf or FALSE to disable the resetting of p-values to NA
  deseq2_results_independentFiltering = TRUE,
  design_var = NULL,
  gene_expression_cols = NULL,
  gmt_file_fdr_cutoffs = c(0.2, 0.05),
  gmt_file_pvalue_cutoffs = c(0.05),
  gmt_ref = "Gene signature created from custom analysis.",
  imported_annotation = NULL,
  my_grouping = NULL,
  output_dir,
  sample_key_col = "Sample_ID"
) {
  
  my_script = "model_differential_expression.R"
  
  my_annotation = paste0(analysis_method," analysis: ", my_script)
  a = function(new_text){
    env = parent.env(environment())
    assign('my_annotation', c(get('my_annotation', envir = env), new_text), envir = env)
  }
  
  
  if(output_dir %>%is.null){
    output_dir = get_analysis_dir("gene_expression_analysis")
    output_dir = file.path(output_dir, analysis_method)
    dir.create(output_dir, showWarnings = FALSE)
  }
  
  if(my_grouping %>% is.null){
    stop("This function requires a grouping by which to split the samples into groups.")
  } 
  
  a(imported_annotation)
  a("")
  
  # dat = import$dat
  # rm(import)
  # cat("Data loading complete.\n")
  
  # if(matched_sample){
  #   # samples that don't have a sample for each level of my_grouping are dropped
  #   n_levels = length(levels(dat[[my_grouping]]))
  #   for(this_individual in levels(dat[[patient_key_col]])){
  #     subdat = dat[dat[[patient_key_col]] == this_individual,  ]
  #     if(nrow(subdat) < n_levels){
  #       dat = dat[dat[[patient_key_col]] != this_individual, ]
  #     }
  #   }
  # }
  
  
  if(base_title %>% is.null){
    base_title = paste0(analysis_method, " Analysis",
                        " by ", my_grouping, '')
  }
  
  if(base_file_name %>% is.null){
    base_file_name = paste0(analysis_method, '_by_', my_grouping)
  }
  
  dir.create(output_dir, showWarnings = FALSE)  # make output folder
  
  gene_expression_cols = gene_expression_cols[gene_expression_cols %ni% sample_key_col]
  
  dat = my_dt %>% as.data.frame(); rm(my_dt)
  
  
  dat = dat[order(dat[,my_grouping]), ]
  
  # drop any extra levels that might be there
  if(class(dat[[my_grouping]]) == "factor"){
    dat[[my_grouping]] = factor(dat[[my_grouping]], levels = levels(dat[[my_grouping]])[levels(dat[[my_grouping]]) %in% dat[[my_grouping]]])
  }
  
  n_classes = length(levels(factor(dat[,my_grouping])))
  
  absent_genes = gene_expression_cols[gene_expression_cols %ni% names(dat)]
  if(length(absent_genes) > 0){
    warning(paste0("These ", length(absent_genes), " genes could not be found in the data:\n", paste0(absent_genes, collapse= ", ")))
    gene_expression_cols = gene_expression_cols[gene_expression_cols %ni% absent_genes]
  }
  
  gene_dat = dat[gene_expression_cols]
  clin_dat = dat[names(dat) %ni% gene_expression_cols]
  
  # DO NOT CHANGE THE ORDER OR NUMBER OF ROWS AFTER THIS POINT!!!!!!!!!!!!!!!!!
  
  display_names = clin_dat$Sample_ID
  
  rownames(gene_dat) = display_names
  
  
  conditions=clin_dat[[my_grouping]]
  
  if(tolower(analysis_method) == "deseq"){
    warning("DESeq has not been debugged yet since major changes were implemented.")
    # a("Processing differential expression using DESeq.")
    # 
    # if(matched_sample){
    #   a(paste0("!!!!! DESeq doesn't support a paired test, so this will proceed with the unpaired analysis.  ",
    #            "Better to use 'SAM' or 'DESeq2' for paired tests. The data will be limited to samples that ", 
    #            "have a match, but the model will not know about their relationship.") %>% wrap_sentence) 
    # }
    # a("FDR corrected pValues are calculated using the Benjamin-Hochberg method.")
    # 
    # countTable = data.frame(t(data.matrix(gene_dat)))
    # rownames(countTable) = colnames(gene_dat)
    # 
    # library(DESeq)
    # cds = newCountDataSet(countTable, conditions)
    # cds = estimateSizeFactors(cds)
    # cds = estimateDispersions(cds)
    # classifier_levels = levels(factor(conditions))
    # res = nbinomTest( cds, classifier_levels[1], classifier_levels[2] )
    # 
    # output_stats = data.frame(
    #   Gene_Name = res$id,
    #   Gene_ID = get_gene_element(res$id, 2),
    #   Fold_Change = res$foldChange,
    #   FDR_pValue = res$padj,
    #   pValue = res$pval
    # )
    # 
    # output_stats = output_stats[ rev(order(output_stats$Fold_Change)), ]
    # 
  } else if (tolower(analysis_method) == "deseq2"){
    cat("Running DESeq2.\n")
    
    # DESeq2 vignette says it's okay to use RSEM input that has been imported through tx import
    #    The tx import DESeq2 approach uses rounded estimated gene counts (but not normalized) 
    #    instead of the raw count of fragments which can be unambiguously assigned to a gene.
    #    This is what we did in assembling the formatted matrix.
    # Important that data are not normalized...
    a("Processing differential expression using DESeq2.")
    a("FDR corrected pValues are calculated using the Benjamin-Hochberg method.")
    
    library(DESeq2)
    countTable = t(data.matrix(gene_dat))
    countTable = apply(countTable, 2, function(x){round(x)})
    
    countTable = data.frame(countTable)
    rownames(countTable) = colnames(gene_dat)
    
    if(!is.null(design_var)){
      stop("design_var has not been setup yet.")
      # Mike Love was the source on doing a paired comparison: https://support.bioconductor.org/p/58893/
      dds = DESeqDataSetFromMatrix(countData = countTable, 
                                   colData = data.frame(conditions, Patient_ID = factor(clin_dat[[patient_key_col]])), 
                                   design = ~ Patient_ID + conditions)
    } else { 
      dds = DESeqDataSetFromMatrix(countData = countTable, 
                                   colData = data.frame(conditions), 
                                   design = ~ conditions)
    }
    
    if (core_number > 1){
      library("BiocParallel")
      register(MulticoreParam(core_number))
    }
    
    if (n_classes == 2){
      my_test = "Wald"
      dds = DESeq(dds, test = my_test, parallel = core_number > 1)
    } else if (n_classes > 2){
      my_test = "LRT"
      dds = DESeq(dds, test = my_test, parallel = core_number > 1, reduced = ~ 1)
    } else {
      stop(paste0("You must have >= 2 groups to run this analysis."))
    }
    
    if ( !is.null(deseq2_results_cooksCutoff)){
      res = results(dds, independentFiltering = deseq2_results_independentFiltering, cooksCutoff = deseq2_results_cooksCutoff)
    } else {
      res = results(dds, independentFiltering = deseq2_results_independentFiltering)
    }
    
    # not sure what to do with LRT data at this point...
    
    
    output_stats = data.frame(
      Gene_Combined_Name = rownames(res),
      Gene_Name = sapply(rownames(res), function(x){strsplit(x, '|', fixed = TRUE)[[1]][1]}) %>% as.character, 
      Gene_ID = sapply(rownames(res), function(x){strsplit(x, '|', fixed = TRUE)[[1]][2]}) %>% as.character,
      Fold_Change = 2^res$log2FoldChange,
      FDR_pValue = p.adjust(res$pvalue, method = "BH"),
      pValue = res$pvalue
    )
    
    # add DeSeq
    deseq2_coef = coef(dds)
    colnames(deseq2_coef) = c("Deseq2_Intercept", "DeSeq2_Coef")
    deseq2_coef = data.frame(Gene_Combined_Name = rownames(deseq2_coef), deseq2_coef)
    rownames(deseq2_coef) = NULL
    output_stats = merge(output_stats, deseq2_coef, by = "Gene_Combined_Name")
    
    
    
    output_stats = output_stats[ rev(order(output_stats$Fold_Change)), ]
    
  } 
  
  # make gmt output files 
  # gmt_file_pvalue_cutoffs = c(0.01),
  # gmt_file_fdr_cutoffs = c(0.2, 0.05)
  cat("Prepping gene signatures.\n")
  id_gmt_file_output = c()
  name_gmt_file_output = c()
  stats_cols = c("pValue", "FDR_pValue")
  for(stats_col in stats_cols){
    if(stats_col == "pValue"){
      gmt_cutoffs = gmt_file_pvalue_cutoffs
    } else if (stats_col == "FDR_pValue"){
      gmt_cutoffs = gmt_file_fdr_cutoffs
    }
    for (gmt_cutoff in gmt_cutoffs){
      sig_stats = output_stats[which(output_stats[,stats_col]<=gmt_cutoff),]
      
      # need to figure out whther these are entrez id's or hgnc
      # we are expecting an hgnc and entrez id seperated by a pipe but need to be ready for anything
      fold_by_col_1_names = sig_stats$Fold_Change
      names(fold_by_col_1_names) = sig_stats$Gene_Name %>% as.character()
      fold_by_col_1_names = fold_by_col_1_names[!is.na(names((fold_by_col_1_names)))]
      fold_by_col_1_names = fold_by_col_1_names[!duplicated(names(fold_by_col_1_names))]
      
      fold_by_col_2_names = sig_stats$Fold_Change
      names(fold_by_col_2_names) = sig_stats$Gene_ID %>% as.character()
      fold_by_col_2_names = fold_by_col_2_names[!is.na(names((fold_by_col_2_names)))]
      fold_by_col_2_names = fold_by_col_2_names[!duplicated(names(fold_by_col_2_names))]
      
      for (fold_by_names in list(fold_by_col_1_names, fold_by_col_2_names)){
        if(sum(!is.na(names(fold_by_names))) > 0){ # if there are any names
          linerized_names = paste0(names(fold_by_names), collapse = "")
          is_id = grepl("^[[:digit:]]+$", linerized_names)
          names_or_ids = ifelse(is_id, "ids", "names")
          gene_set_base_name = paste0(base_file_name, "__", stats_col, "_", gmt_cutoff)# , "_sig_genes"
          up_genes = names(fold_by_names[fold_by_names > 1])
          down_genes = names(fold_by_names[fold_by_names < 1])
          if(length(up_genes) > 0){
            gene_set_name = paste0(gene_set_base_name, "__up")
            if(is_id){
              id_gmt_file_output = c(id_gmt_file_output, paste0(c(gene_set_name, gmt_ref, paste0(up_genes, collapse = "\t")), collapse = "\t"))
            } else {
              name_gmt_file_output = c(name_gmt_file_output, paste0(c(gene_set_name, gmt_ref, paste0(up_genes, collapse = "\t")), collapse = "\t"))
            }
          }
          if(length(down_genes) > 0){
            gene_set_name = paste0(gene_set_base_name, "__down")
            if(is_id){
              id_gmt_file_output = c(id_gmt_file_output, paste0(c(gene_set_name, gmt_ref, paste0(down_genes, collapse = "\t")), collapse = "\t"))
            } else {
              name_gmt_file_output = c(name_gmt_file_output, paste0(c(gene_set_name, gmt_ref, paste0(down_genes, collapse = "\t")), collapse = "\t"))
            }
          }
        }
      }
    }
  }
  if(length(id_gmt_file_output) > 0){
    gmt_id_path = file.path(output_dir, paste0(base_file_name, "_gene_signature_ids.gmt.txt"))
    writeLines(id_gmt_file_output, gmt_id_path)
    cat("Producing a gmt file with gene ids.\n")
  } else {
    gmt_id_path = NULL
    cat("No significant genes were found to make a gene signature with gene ids. Either no genes were significant or only gene names were provided.\n")
  }
  
  if(length(name_gmt_file_output) > 0){
    gmt_symbol_path = file.path(output_dir, paste0(base_file_name, "_gene_signature_symbols.gmt.txt"))
    writeLines(name_gmt_file_output, gmt_symbol_path)
    cat("Producing a gmt file with gene names.\n")
  } else {
    gmt_symbol_path= NULL
    cat("No significant genes were found to make a gene signature with gene names. Either no genes were significant or only gene ids were provided.\n")
  }
  
  
  a("Output stats with a pvalue <= 0.05."); a("")
  names(output_stats) %<>% gsub(".", "_",., fixed = T)
  
  #sig_stats = output_stats[which(output_stats$pValue<=0.05),] # for writing to text file
  # rownames(output_stats) = output_stats$Gene_Name
  # heatmap_pValues = output_stats[names(gene_dat), "pValue"]
  # names(heatmap_pValues) = names(gene_dat)
  # heatmap_FDR_pValues = output_stats[names(gene_dat), "FDR_pValue"]
  # names(heatmap_FDR_pValues) = names(gene_dat)
  # heatmap_fold_change = output_stats[names(gene_dat), "Fold_Change"]
  # names(heatmap_fold_change) = names(gene_dat)
  # 
  # sig_genes = sig_stats$Gene_Name # get the list of genes that were significant in order of high to low fold change
  # #     heatmap_stats = data.frame(feature = pretty_gene_name(sig_stats$Gene_ID %>% as.character), pValue = sig_stats$FDR_pValue/100)
  # #     gene_dat = gene_dat[ , sig_genes %>% as.character] 
  # 
  # sig_stats$Gene_ID = sapply(sig_stats$Gene_Name %>% as.character, function(x){strsplit(x, '|', fixed = TRUE)[[1]][2]}) %>% as.character
  # sig_stats$Gene_Name = sapply(sig_stats$Gene_Name %>% as.character, function(x){strsplit(x, '|', fixed = TRUE)[[1]][1]}) %>% as.character
  # 
  # # output stats of sam analysis
  stats_path = file.path(output_dir, paste0(base_file_name, ".stats.tsv"))
  fwrite(output_stats[order(output_stats$pValue), ], stats_path, sep = "\t")
  
  # heatmap_pack = list(sample_names = clin_dat[[sample_key_col]],
  #                     classifiers = conditions,
  #                     input_df = gene_dat, # +1 to make sure 0's aren't turned to NA's
  #                     base_title = base_title,
  #                     base_file_name = base_file_name,
  #                     output_dir = output_dir,
  #                     # pdf_metadata = pdf_metadata,
  #                     # argument_list = argument_list,
  #                     my_model = analysis_method,
  #                     annotation = c(my_annotation),
  #                     FDR_pValue = heatmap_FDR_pValues,
  #                     pValue = heatmap_pValues,
  #                     fold_change = heatmap_fold_change
  # )
  
  if(n_classes == 2){ 
    
    volcano_pack = list(
      base_file_name = base_file_name,
      # base_title = base_title,
      colname_of_feature_groups = NULL,
      colname_of_feature_names = "Gene_Name",
      colname_of_fold_change = "Fold_Change",
      colname_of_pValue = "pValue",
      colname_of_FDR = "FDR_pValue",
      my_annotation = my_annotation,
      my_dt = output_stats %>% as.data.table(),
      ordered_factors = levels(conditions)[1:2],
      output_dir = output_dir
    )
    
  } else {
    a("There were more than 2 categories so a volcano plot wasn't possible. Compare two categories at a time for this." %>% wrap_sentence)
    volcano_pack = NULL 
  }
  
  a(paste0(" End of ", analysis_method," analysis") %>% as.footer)
  
  annotation_path = file.path(output_dir, paste0(base_file_name, get_note_extension()))
  write.table(my_annotation, annotation_path, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  
  # rerun_path = file.path(output_dir, paste0(base_file_name, get_config_extension()))
  # write.table(config_script, rerun_path, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  
  return(list(volcano_pack = volcano_pack, gmt_id_path = gmt_id_path, gmt_symbol_path = gmt_symbol_path, stats_path = stats_path))
  
  
} # end gene_expression_analysis
