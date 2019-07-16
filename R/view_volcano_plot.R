# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# view_volcano_plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#
#
#
#' @title Create a volcano plot of analysis results.
#'
#' @description
#' Fold change and significance are used to produce a volcano plot of data for gene expression, 
#' DAVID analysis and gene set colloctions
#' 
#' @details
#' \code{bgv_volcano_plot} points may be labeled based on being over a certain significance or just
#' by labeling a certain number to label.  Clusters of points may also be colored to show how the group is doing.
#'
#'   
#' @param my_pack A named list for analysis modules to send parameters to the view function in a 
#' more compact fashion:
#' \itemize{
#'   \item base_file_name
#'   \item base_title
#'   \item colname_of_feature_groups
#'   \item colname_of_feature_names
#'   \item colname_of_fold_change
#'   \item colname_of_pValue
#'   \item colname_of_FDR
#'   \item imported_annotation
#'   \item my_dt
#'   \item ordered_factors
#'   \item output_dir
#' }
#' These list items may be overridden by specifying them individually.  They are described in 
#' more detail below.
#' 
#' @param base_file_name Overrides \code{my_pack$base_file_name}. String to specify the file name.
#' @param base_title Overrides \code{my_pack$base_title}. String to specify the title of the graph. 
#' @param colname_of_feature_groups Overrides \code{my_pack$colname_of_feature_groups}. String indicating the name of the column to get the feature groups from.  
#'   If providing groups the points will be colored according to these groups and not significance.
#' @param colname_of_feature_names Overrides \code{my_pack$colname_of_feature_names}. String indicating the name of the column to pull the feature names from.
#' @param colname_of_fold_change Overrides \code{my_pack$colname_of_fold_change}. String indicating the name of the column to get the fold change from.
#' @param colname_of_pValue Overrides \code{my_pack$colname_of_pValue}. String indicating the name of the column to get the pValue from.
#' @param colname_of_FDR Overrides \code{my_pack$colname_of_FDR}. String indicating the name of the column to get the FDR from. If not provided the fdr_method 
#'   will be used to calculate FDR values.  fdr_method should be set to "none" if no correction is desired.
#' @param imported_annotation Overrides \code{my_pack$my_annotation}. Annotations describe the 
#'   steps taken during the analysis module and previous steps. This function will also add 
#'   annotation steps to provide a full description of how the graph was generated.
#' @param my_dt Data table to provide the information for the graph.  It should have columns for pValue and fold_change 
#'   as a minimum.
#' @param ordered_factors Overrides \code{my_pack$ordered_factors}. A character string of length 2 to tell the x axis what it should be labeled
#' @param output_dir Overrides \code{my_pack$output_dir}. Path to folder where the heatmap output should go.  
#'   Analysis modules will determine this automatically. 
#'   
#' @param add_pdf_metadata Boolean indicating whether the pdf metadata should be chnaged.  Is pdftk present?
#' @param axis.title.x.size Numeric indicating the size of the x axis title.
#' @param axis.title.y.size Numeric indicating the size of the y axis title.
#' @param FDR_label_threshhold Number to specify threshold over which points will not be labeled.
#' @param fdr_method String indicting the type of fdr method to use if no colname_of_FDR is given.  See \code{\link{p.adjust}} 'method'. 
#' @param FDR_sig_alpha numeric to indicate the alpha value of fdr significant points
#' @param label_font_size Numeric to indicate the size of the point labels.
#' @param label_pValue_line Boolean to determine whether the pValue line should be labeled.
#' @param log_completion Binary indicting whether the graph should be added to the table of completed analyses.
#' @param max_number_of_labels integer indicating the max number of points that will be labeled.
#' @param nonsig_label String to label legend for non-significant points.
#' @param num_top_groups_to_color Integer that tells how many of the top groups (using group argument) should be labeled and colored.
#' @param n_to_label Integer to specify how many points should be labeled.  More than 25 or so can get a bit cluttered. This argument supercedes FDR and pValue label threshholds. Set to 0 if you don't want labels
#' @param output_subfolder String that specifies what folder the output will go into in the main output directory. 
#' @param pdf_bookmark_title String indicating the name of the bookmark that will be used when concatentaing the pdf's into a report
#' @param pdf_height Sets the height (inches) of the pdf output
#' @param pdf_width Sets the width (inches) of the pdf output
#' @param plot_by_FDR boolean to indicate if the y axis should be against the FDR_corrected pvalues
#' @param pValue_label_threshhold Number to specify threshold over which points will not be labeled.
#' @param pValue_line Number to specify where to draw a line to label a specific pValue. Usually set to 0.05 if used.
#' @param sig1_cutoff Numeric to indicate the upper bound for the most significant points.
#' @param sig1_label String to label legend for the most significant points.
#' @param sig1_type Either "FDR" or "pValue". Indicates what should be used to determine the most significant points.
#' @param sig2_cutoff Numeric to indicate the upper bound for the second most significant points.
#' @param sig2_label String to label legend for the second most significant points.
#' @param sig2_type Either "FDR" or "pValue". Indicates what should be used to determine the second-most significant points.
#' @param strings_to_cut_from_labels In certain gene group names a lot of the names can be very long and not very informative. This vector of character strings allows you to specify strings that will be cut out of the names. Leave off leading and trailing "_".
#' @param x_axis_label String to specify label for x axis.
#' 
#' @return
#'   1) Volcano plot; 
#'   2) Volcano plot pdf;
#'     
#' @section Todos:
#' \itemize{
#'   \item Debug for groups
#' }
#' 
#' @section Future changes:
#' \itemize{
#'   \item ...
#' }
#' 
#' @family view
#' 
#' @export
view_volcano_plot = function(# passed from model app
  my_pack = NULL,
  base_file_name = NULL,
  base_title = NULL,
  colname_of_feature_groups = NULL,
  colname_of_feature_names = NULL,
  colname_of_fold_change = NULL,
  colname_of_pValue = NULL,
  colname_of_FDR = NULL,
  imported_annotation = NULL,
  my_dt = NULL,
  ordered_factors = NULL,
  output_dir = NULL,
  # ----------
  add_pdf_metadata = TRUE,
  axis.title.x.size = 12,
  axis.title.y.size = 16,
  FDR_label_threshhold = NULL,
  fdr_method = "BH",
  FDR_sig_alpha = 1,
  label_font_size = 2.5,
  label_pValue_line = TRUE,
  log_completion = FALSE,
  max_number_of_labels = 50,
  nonsig_label = NULL,
  num_top_groups_to_color = 6,
  n_to_label = NULL,
  output_subfolder = "XX_Volcano_Plot",
  pdf_bookmark_title = "02_Volcano_Plot",
  pdf_height = 7,
  pdf_width = 9,
  plot_by_FDR = FALSE,
  pValue_label_threshhold = NULL,
  pValue_line = NULL,
  sig1_cutoff = 0.05,
  sig1_label = NULL,
  sig1_type = "FDR",
  sig2_cutoff = 0.05,
  sig2_label = NULL,
  sig2_type = "pValue",
  strings_to_cut_from_labels = c("MOUSECYC_MM", "PWY", "WIKIPATHWAYS_MM", "PANTHER_MM", "BIOCARTA_MM", 
                                 "WIKIPATHWAYS_MM", "KEGG_MM", "TF_MM", "REACTOME_MM", "TFACTS_MM", 
                                 "PID_", "KEGG_", "BIOCARTA", "REACTOME_", "NABA_"),
  x_axis_label = NULL
){
  
  
  my_script = "view_volcano_plot.R"
  my_annotation = paste0("Produce volcano plot: ", my_script)
  a = function(new_text){
    env = parent.env(environment())
    assign('my_annotation', c(get('my_annotation', envir = env), new_text), envir = env)
  }
  
  if(my_pack %>% is_not_null){
    if(base_file_name %>% is.null){ base_file_name = my_pack$base_file_name }
    if(base_title %>% is.null()){ base_title = my_pack$base_title }
    if(colname_of_feature_groups %>% is.null){ colname_of_feature_groups = my_pack$colname_of_feature_groups }
    if(colname_of_feature_names %>% is.null){ colname_of_feature_names = my_pack$colname_of_feature_names }
    if(colname_of_fold_change %>% is.null){ colname_of_fold_change = my_pack$colname_of_fold_change }
    if(colname_of_pValue %>% is.null){ colname_of_pValue = my_pack$colname_of_pValue }
    if(colname_of_FDR %>% is.null){ colname_of_FDR = my_pack$colname_of_FDR }
    if(imported_annotation %>% is.null){ imported_annotation = my_pack$my_annotation }
    if(my_dt %>% is.null){ my_dt = my_pack$my_dt }
    if(ordered_factors %>% is.null){ ordered_factors = my_pack$ordered_factors }
    if(output_dir %>% is.null){ output_dir = my_pack$output_dir }
  }
  
  # let's do this as a data.frame for now
  input_df = my_dt %>% as.data.frame(); rm(my_dt)
  
  if(is.null(pdf_bookmark_title) || is.na(pdf_bookmark_title)){
    pdf_bookmark_title = gsub("_", " ", base_file_name)
  }
  
  if(!is.null(n_to_label) && is.na(n_to_label)) n_to_label = NULL # double && important here as null will crash is.na
  if(!is.null(pValue_label_threshhold) && is.na(pValue_label_threshhold)) pValue_label_threshhold = NULL
  if(!is.null(FDR_label_threshhold) && is.na(FDR_label_threshhold)) FDR_label_threshhold = NULL
  
  # improt annotation if there
  if(!is.null(imported_annotation)){ a(format_imported_annotation(imported_annotation)) }

  # add subfolder if there
  if(output_subfolder %>% is_not_null){ output_dir = file.path(output_dir, output_subfolder) }
  dir.create(output_dir, showWarnings = FALSE)  # make output folder
  
  
  
  # need error if any of the critical columns are NULL
  required_columns = c("colname_of_fold_change", "colname_of_pValue")
  for (required_column in required_columns){
    my_col = eval(parse(text = required_column))
    if(is.null(my_col)){
      stop(paste0(my_col, "was NULL.  This is a required column."))
    }
    if(my_col %ni% names(input_df)){
      stop(paste0(required_column, " ('",my_col, "') was not found in my_dt. This column is required."))
    } else {
      # drop features if they are missing critical values
      missing_data = !complete.cases(input_df[[my_col]])
      if(sum(missing_data) > 0){
        my_warning = paste0(sum(missing_data), " were missing data for ", my_col, " and were dropped.")
        warning(my_warning)
        a(my_warning)
        input_df = input_df[!missing_data, ]
      }
    }
  }

  # now need error if any non critical column names are missing
  other_columns = c("colname_of_feature_groups", "colname_of_feature_names", "colname_of_FDR")
  for (other_column in other_columns){
    my_col = eval(parse(text = other_column))
    
    if(!is.null(my_col) && my_col %ni% names(input_df)){
      stop(paste0(other_column, " ('",my_col, "') was not found in my_dt. Did you specify the correct name?"))
    }
  }
  
  
  # standardize these names
  volcano_df = data.frame(
    Feature_Name = input_df[[colname_of_feature_names]], 
    Fold_Change = input_df[[colname_of_fold_change]], 
    pValue = input_df[[colname_of_pValue]]
  )
  
  
  if(is.null(colname_of_FDR)){
    volcano_df$FDR = p.adjust(input_df[[colname_of_pValue]], method = fdr_method)
  } else {
    volcano_df$FDR = input_df[[colname_of_FDR]]
  }
  
  
  if(!is.null(colname_of_feature_groups)){
    volcano_df$Group = input_df[[colname_of_feature_groups]]
  }
  
  rm(input_df)
  
  
  #output_df = volcano_df
  max_significance = pvalue_stars(min(volcano_df$FDR, na.rm = T))
  
  # make legend labels
  if(sig1_label %>% is.null()){
    sig1_label = paste0(sig1_type, " <= ", sig1_cutoff)
  }
  if(sig2_label %>% is.null()){
    sig2_label = paste0(sig2_type, " <= ", sig2_cutoff)
  }
  if(nonsig_label %>% is.null()){
    nonsig_label = paste0(sig2_type, " > ", sig2_cutoff)
  }


  
  FDR_string = sig1_label
  unadjusted_string = sig2_label
  ns_string = nonsig_label
  
  if(sig1_type %ni% c("FDR", "pValue")){
    warning(paste0("sig1_type, ", sig1_type, ", is not recognized.  It should be either 'FDR' or 'pValue'."))
  }
  
  if(sig2_type %ni% c("FDR", "pValue")){
    warning(paste0("sig2_type, ", sig2_type, ", is not recognized.  It should be either 'FDR' or 'pValue'."))
  }
  
  volcano_df$Significance = ifelse(volcano_df[[sig1_type]]<=sig1_cutoff, 
                                   FDR_string,
                                   ifelse(volcano_df[[sig2_type]]<=sig2_cutoff,
                                          unadjusted_string,
                                          ns_string))
  sig_levels = c(FDR_string, unadjusted_string, ns_string)
  sig_colors = c( "red", "darkred", "black")
  volcano_df$Significance = factor(volcano_df$Significance, levels = sig_levels)
  

  if (plot_by_FDR){
    y_axis_label = "-Log10 FDR-Corrected pValue"
    plot_column = "FDR"
  } else {
    y_axis_label = "-Log10 Unadjusted pValue"
    plot_column = "pValue"
  }
  
  
  if(x_axis_label %>% is.null() && length(ordered_factors) == 2){
    x_axis_label = paste0("<-- ", ordered_factors[1], " high     |     Log2(Fold Change)     |     ", ordered_factors[2], " high -->")
  } else { 
    x_axis_label = paste0("Log2(Fold Change)")
  }
  
  if(colname_of_feature_groups %>% is.null){
    
    color_values = character(0)
    for (sig_index in 1:length(sig_levels)) {
      this_level = sig_levels[sig_index]
      if (this_level %in% volcano_df$Significance){
        color_values = c(color_values, sig_colors[sig_index])
      }
    }
    
    
    volcano_plot = ggplot(data = volcano_df, aes(x = log2(Fold_Change), y = -log10(eval(parse(text=plot_column))), col = Significance, 
                                                 alpha = ifelse(Significance == FDR_string, FDR_sig_alpha, 0.1))) + 
      scale_alpha(range=c(0.35,1), limits=c(0,1)) + guides(alpha=FALSE) +
      geom_point() + scale_colour_manual(name="Sig.",  values =color_values) +
      scale_y_continuous(expand = c(0.01, 0))
    
    
  } else {
    
    volcano_df$Feature_Name = paste0("Group ", volcano_df$Group, ": ", volcano_df$Feature_Name)
    volcano_df$Group = volcano_df$Group %>% as.character %>% as.numeric
    volcano_df$Group[ as.numeric(as.character(volcano_df$Group)) > num_top_groups_to_color] = NA
    volcano_df$Group = factor(volcano_df$Group)
    
    volcano_plot = ggplot(data = volcano_df, aes(x = log2(Fold_Change), y = -log10(eval(parse(text=plot_column))), color = Group,  
                                                 alpha = ifelse(Significance == FDR_string, FDR_sig_alpha, 0.1))) + 
      geom_point() +
      scale_alpha(range=c(0.35,1), limits=c(0,1)) + guides(alpha=FALSE, color = FALSE) + 
      scale_colour_hue(na.value = "black")
    
  }
  
  # add labels
  volcano_plot = volcano_plot + labs(title = paste0(base_title), x = x_axis_label, y = y_axis_label) + 
                          theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size=22, margin=margin(0,0,20,0)),
                          axis.title.x = element_text(size = axis.title.x.size, margin=margin(15,0,0,0)),
                          axis.title.y = element_text(size = axis.title.y.size, margin=margin(0,15,0,0)),
                          axis.text = element_text(size = 12),
                          legend.text = element_text(size = 12),
                          legend.title = element_text(size = 18)
  ) 
  
  
  if(pValue_line %>% is_not_null){
    #http://stackoverflow.com/questions/7705345/how-can-i-extract-plot-axes-ranges-for-a-ggplot2-object
    #ggplot version was changed by bioconductor changed and it affected this 
    # x_range = ggplot_build(volcano_plot)$panel$ranges[[1]]$x.range
    # y_range = ggplot_build(volcano_plot)$panel$ranges[[1]]$y.range    
    x_range = ggplot_build(volcano_plot)$layout$panel_ranges[[1]]$x.range
    y_range = ggplot_build(volcano_plot)$layout$panel_ranges[[1]]$y.range
    
    volcano_plot = volcano_plot + geom_hline(yintercept = -log10(pValue_line), linetype="dotted")
    if (label_pValue_line){
      volcano_plot = volcano_plot + annotate("text", x = x_range[1], y = -log10(pValue_line), 
                                             label = paste0("pValue = ", pValue_line), hjust = 0.3, vjust = 1.2)
    }
  }

  labeled_points = volcano_df

  if(n_to_label %>% is_not_null){
    if( n_to_label > 0){
      if(n_to_label > nrow(labeled_points)){
        n_to_label = nrow(labeled_points)
        a(""); a(paste0("Labeled all points."))
      } else {
        a(""); a(paste0("Labeled top ", n_to_label," points."))
      }
      labeled_points = labeled_points[order(labeled_points$pValue)[1:n_to_label], ]
    } else {
      a(""); a(paste0("No points labeled."))
      labeled_points = labeled_points[0, ]
    }
  } else {
    
    
    if(pValue_label_threshhold %>% is.null() & FDR_label_threshhold %>% is.null() ){
      if(sum(labeled_points$FDR <=  0.05, na.rm = T) > 0){
        FDR_label_threshhold = 0.05
      } else if(sum(labeled_points$pValue <=  0.05, na.rm = T) > 0){
        pValue_label_threshhold = 0.05
      }
    }
    
    if(pValue_label_threshhold %>% is_not_null){
      a("")
      a(paste0("Labeled items with pValue <= ", pValue_label_threshhold, "."))
      labeled_points = labeled_points[labeled_points$pValue <= pValue_label_threshhold, ]
    } else if(FDR_label_threshhold %>% is_not_null){
      a("")
      a(paste0("Labeled items with FDR corrected pValue <= ", FDR_label_threshhold, "."))
      labeled_points = labeled_points[labeled_points$FDR <= FDR_label_threshhold, ]
    } 
    
    if(nrow(labeled_points) > max_number_of_labels){
      labeled_points = labeled_points[order(labeled_points$pValue)[1:max_number_of_labels], ]
    }
  }
  

  if(nrow(labeled_points) > 0){
    library(ggrepel)
    
    # trim text
    if(strings_to_cut_from_labels %>% is_not_null){
      for(this_str in strings_to_cut_from_labels){
        labeled_points$Feature_Name = gsub(this_str, "", labeled_points$Feature_Name)
        labeled_points$Feature_Name = gsub("__", "_", labeled_points$Feature_Name) # get rid of doubles created by removals
      }
      
      labeled_points$Feature_Name = gsub("^_", "", labeled_points$Feature_Name) # get rid of doubles created by removals
      labeled_points$Feature_Name = gsub("_$", "", labeled_points$Feature_Name) # get rid of doubles created by removals
      
    }
    

    if(colname_of_feature_groups %>% is.null){
       volcano_plot = volcano_plot + geom_label_repel(data = labeled_points, aes(label=Feature_Name), alpha = 1, 
                                                       color = 'black', size = label_font_size, show.legend = FALSE, 
                                                       # seems like color of font and outline are linked
                                                       box.padding = unit(0.3, "lines"),# sets how spread apart the boxes are
                                                       point.padding = unit(0.3, "lines"),# sets how far the lines sits away from the point
                                                       label.padding = unit(0.15, "lines"), # sets how big the box around the text is
                                                       segment.color = "#666666", segment.size = 0.25,
                                                       force = 2)
    } else { # for DAVID clustered analysis
      volcano_df$Significance[volcano_df$Group %>% is.na] = ns_string
      
      volcano_plot = volcano_plot + geom_label_repel(data = labeled_points, aes(label=Name), alpha = 1, 
                                                     size = label_font_size, show.legend = FALSE, 
                                                     # seems like color of font and outline are linked
                                                     box.padding = unit(0.3, "lines"),# sets how spread apart the boxes are
                                                     point.padding = unit(0.3, "lines"),# sets how far the lines sits away from the point
                                                     label.padding = unit(0.15, "lines"), # sets how big the box around the text is
                                                     segment.color = "#666666", segment.size = 0.25,
                                                     force = 2)
    }
  }
  volcano_plot_gt <- ggplot_gtable(ggplot_build(volcano_plot))
  volcano_plot_gt$layout$clip[volcano_plot_gt$layout$name=="panel"] <- "off"
  
  
  save_path = file.path(output_dir, paste0(base_file_name,".pdf"))
  ggsave(save_path, volcano_plot_gt, device = "pdf", width = pdf_width, height = pdf_height)
  
  if(add_pdf_metadata){
    # pdf_metadata is still being passed and could be added if we want, but for now just adding what we need for latex to make the bookmarks
    # format is a string c("value_id1: value1", "value_id=2: value2" 
    # to combine into nice pdfs we need the pdf_width, pdf_height, significance,	pdf_bookmark_title	pdf_bookmark_id
    
    # for significance if an fdr correction is shown we will output the highest number of stars
    
    my_meta_data = c(
      paste0("pdf_bookmark_title: ", pdf_bookmark_title),
      paste0("significance: ", max_significance)
    )
    
    write_metadata(save_path, my_meta_data)
    
  }
  
  
  #write_metadata(save_path, pdf_metadata)
  if(log_completion){
    warning("Haven't fixed toolkit's add_to_completed_analyses yet")
    
#     add_to_completed_analyses(
#       classifier = models_argument_list$classifier,
#       included = models_argument_list$inclusion_list,
#       excluded = models_argument_list$exclusion_list,
#       model = my_model,
#       view = "Volcano_Plot",
#       path = save_path)
  }
  a(' Volcano plot done' %>% as.footer)
  annotation_path = file.path(output_dir, paste0(base_file_name, get_note_extension()))
  write.table(my_annotation, annotation_path, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  
  # save the data that is on the plot
  #data_shown_df = data.frame(Gene = name)
  
  
  return(volcano_plot)
  
}
