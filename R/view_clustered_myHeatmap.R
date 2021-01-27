# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# view_clustered_myHeatmap
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#
#
#' @title Create a heatmap of analysis results.
#'
#' @description
#' Produces a heatmap of scaled values by row and then makes an unscaled average value column to the right of it.
#' 
#' @details
#' \code{view_clustered_myHeatmap} was created to make heatmaps of gene expression data, although it could
#' be used to display the distribution of any numeric features.
#' It's important to not scale data prior to heatmap to not disrupt average values.
#'
#' @param my_pack A named list for analysis modules to send parameters to the view function in a 
#' more compact fashion:
#' \itemize{
#'   \item sample_names
#'   \item groups
#'   \item input_df
#'   \item base_file_name
#'   \item base_title
#'   \item output_dir
#'   \item pdf_metadata
#'   \item argument_list
#'   \item my_model
#'   \item annotation
#'   \item FDR_pValue
#'   \item pValue
#'   \item fold_change
#' }
#' These list items may be overridden by specifying them individually.  They are described in 
#' more detail below.
#' 
#' @param sample_names Overrides \code{my_pack$sample_names}. Should be the sample names you want to have for 
#'   the scaled values. They should be in the order you want them displayed. They should correspond to the rows 
#'   of input_df. They should have the same order as thier groups.
#' @param groups Overrides \code{my_pack$groups} A factored vector indicating what groups the 
#'   samples belong to.  Should be the same length as \code{sample_names}. Their order should match sample_names.
#'   They should be factored with the order matching the order in which they should be displayed. They should be 
#'   formatted how you would like them displayed.
#' @param input_df Overrides \code{my_pack$input_df} A data frame with one row per sample.  Column names should
#'   correspond to the features along the y axis of the heatmap. Row order should correspond to the samples/groups
#'   should be in the order in which you would like them displayed.
#' @param base_file_name Overrides \code{my_pack$base_file_name}. String to specify the file name.
#' @param base_title Overrides \code{my_pack$base_title}. String to specify the title of the graph.
#' @param output_dir Overrides \code{my_pack$output_dir}. Path to folder where the heatmap output should go.  
#'   Analysis modules will determine this automatically. 
#' @param pdf_metadata Overrides \code{my_pack$pdf_metadata}. Contains data from the arguments list that should
#'   be put in the metadata of the pdf.
#' @param models_argument_list Overrides \code{my_pack$argument_list}.  Lists all of the arguments sent to the 
#'   analysis module.
#' @param my_model Overrides \code{my_pack$my_model}. Links to all graphs are logged on a d3.js table. This variable 
#'   will specify what model name was used to generate the graph on that table.
#' @param imported_annotation Overrides \code{my_pack$annotation}. Annotations describe the 
#'   steps taken during the analysis module and previous steps. This function will also add 
#'   annotation steps to provide a full description of how the graph was generated.
#' @param FDR_pValue Overrides \code{my_pack$FDR_pValue}. Numeric vector with length equal to \code{ncol(input_df)}. 
#'   If \code{FDR_threshhold} is not \code{NULL} then only features with an \code{FDR_pValue} less than the 
#'   \code{FDR_threshhold} will be included.
#' @param pValue Overrides \code{my_pack$pValue}. Numeric vector with length equal to \code{ncol(input_df)}.
#'   If \code{pValue_threshhold} is not \code{NULL} then only features with a \code{pValue} less than the 
#'   \code{pValue_threshhold} will be included.
#' @param fold_change Overrides \code{my_pack$fold_change}. Numeric vector with length equal to \code{ncol(input_df)}.
#'   If \code{order_by_fold_change == TRUE} then features will be ordered according to this value from high to low, where
#'   'higher' fold change indicates group2 was higher than group1.
#'   
#' @param order_by_fold_change Logical indicating whether the features should be ordered by fold change. See fold_change argument.
#' @param base_for_log_transform Integer indicating the base to transform the data, which will be shown in the averages. 
#'   \code{NULL} will results in no transformation. 
#' @param sample_panel_spacing \code{unit()} used to set the spacing between individual sample panels.
#' @param average_panel_spacing \code{unit()} used to set the spacing between average sample panels, usually 0.
#' @param show_sample_names Boolean indicating whether individual sample names should be listed along the x axis. 
#'   This axis can get cluttered quickly if there are a lot of samples.
#' @param append_classifier_to_sample_name Boolean indicating whether the sample group should be added to the end of the sample names.  
#'   This option isn't really needed after having added facet titles.
#' @param panel_label_switch 'x' switches the x axis label to the bottom. \code{NULL} leaves the labels where they are.
#' @param panel.strip.background For specifying the element_rect of the panel strip background (eg. \code{element_rect(fill = 'gray95')}). 
#' @param panel.strip.text.x For specifying the element_rect of the panel strip (eg. \code{element_text(size=18, angle = 0, hjust = 0.5, vjust = 0)}.
#' @param FDR_threshhold If not \code{NULL}, features with an \code{FDR_pValue} above this value will not be shown.
#' @param pValue_threshhold If not \code{NULL}, features with an \code{pValue} above this value will not be shown.
#' @param display_all_names Boolean indicating whether the feature names should be displayed along the y axis.
#' @param graph_rel_widths Numeric vector of length 2 (eg. \code{c(8, 2.5)}) whose ratio determines the size of the individual sample panel vs the 
#'   averages.
#' @param axis.text.x element_text() object to define the x axis for the sample and average plots
#' @param axis.text.y.size Number indicating the size of the text for the feature names along the y axis.
#'    
#' @param pdf_width Sets the width (inches) of the pdf output
#' @param pdf_height Sets the height (inches) of the pdf output
#' @param log_completion Binary indicting whether the graph should be added to the table of completed analyses.
#' @param primary_title_options List of options to modify primary title. Title content should be set with base_title. eg: list(x=0.5,y=0.965,size=8)
#' @param secondary_title List of options to modify and set secondary title.  eg: list(x=1,y=0.99,size=3,label="FDR < 0.05 shown")
#' @param plot_scale numeric to change how much of the page the plots take up
#' @param legend_title \code{element_text()} object to modify the disply of the legend titles
#' @param legend_text \code{element_text()} object to modify the disply of the legend text
#' @param legend_key_width \code{unit()} for the width of both legends
#' @param legend_key_height \code{unit()} for the height of both legends
#' @param ave_legend_margin \code{ margin(t, r, b, l, "pt")} 
#' @param zscore_legend_margin \code{ margin(t, r, b, l, "pt")}
#' @param ave_plot_margin \code{unit(t, r, b, l, "pt")}
#' @param zscore_plot_margin \code{unit(t, r, b, l, "pt")}
#' @param left_legend_pos string for position of zscore plot
#' @param right_legend_pos string for position of average plot
#' @param axis_title_x \code{element_text()} or \code{element_blank()} for axis of zscore plot
#' @param axis_title_y \code{element_text()} or \code{element_blank()} for axis of zscore plot
#' @param axis_title_x_lab string for axis of zscore plot, usually blank
#' @param axis_title_y_lab string for aixs of zscore plot, usually blank
#' @param zscore_colors vector of colors for z scale, overrides average colors if \code{average_scaled_vlaues == T}
#' @param average_colors vecotr of colors
#' @param average_scaled_values boolean indicating whether the average should be done on the raw gene expression or scaled values
#' @param axis_text_y_remove_strings Strings in this character vector will be removed from the y axis text.
#' @param axis_text_y_nchar Integer to specify the number of characters that will be displayed in the y axis text before clipping
#' @param pdf_bookmark_title String indicating the name of the bookmark that will be used when concatentaing the pdf's into a report
#' @param annotate_pdf Boolean indicating whether the pdf metadata should be chnaged.  Is pdftk present?
#' @param output_subfolder Path to the output sub directory folder. Multiple views of the same analysis should have unique names lest 
#'   they overwrite one another. XX indicates this folder will be dropped from the TOC and bookmarks
#' @param tile_width_factor Changes the size of the tiles to get rid of gaps or create them. Goes into geom_tile width parameter.
#' @param tile_height_factor Changes the size of the tiles to get rid of gaps or create them. Goes into geom_tile height parameter.
#' 
#' @return Several outputs are produced: 
#'   1) Heatmaps; 
#'   2) Heatmaps pdf;
#'   3) An annotation file describing the steps taken on the data from start to finish;
#'   4) A matrix file containing the sample z-scores for each sample and the average values displayed in the heatmap;
#'   5) If \code{log_completion == TRUE}, a link is made to the pdf on the d3.js html indexing 
#'     table using \code{add_to_completed_analyses}.
#'     
#' @section Limitations:
#' \itemize{
#'   \item ...
#' }
#' 
#' @section Future changes:
#' \itemize{
#'   \item ...
#' }
#' 
#' @export
view_clustered_myHeatmap = function(# passed from analysis app
  my_pack = NULL,
    # for overriding my_pack contents------
    sample_names = NULL,
    groups = NULL,
    input_df = NULL,
    base_title = NULL,
    base_file_name = NULL,
    output_dir = NULL,
    pdf_metadata = NULL,
    models_argument_list = NULL,
    imported_annotation = NULL,

  # order_by_fold_change = TRUE,
  # base_for_log_transform = NULL,
  # sample_panel_spacing = unit(0.05, "inches"),
  # average_panel_spacing = unit(0.0, "inches"),
  # show_sample_names = FALSE,
  # append_classifier_to_sample_name = FALSE,
  # panel_label_switch = NULL,
  # panel.strip.background = element_rect(fill = 'gray95'),
  # panel.strip.text.x = element_text(size=18, angle = 0, hjust = 0.5, vjust = 0),
  # FDR_threshhold = NULL,
  # pValue_threshhold = 0.05,
  # display_all_names = TRUE,
  # graph_rel_widths = c(8, 2.5),
  # axis.text.x = element_text(size=12, angle = 90, hjust = 0, vjust = 0.5),
  # axis.text.y.size = NULL,
  # pdf_width = 8,
  # pdf_height = 10,
  # log_completion = FALSE,
  # primary_title_options = list(x=0.5,y=0.965,size=8),
  # secondary_title = list(x=1,y=0.98,size=3, label=""),
  # plot_scale = 1,
  # legend_title = element_text(size=12, hjust = 0, vjust = 0.5),
  # legend_text = element_text(size=10, hjust = 0, vjust = 0.5),
  # legend_key_width = unit(0.1, "inches"),
  # legend_key_height = unit(2, "inches"),
  # ave_legend_margin = margin(0, 0, 0, 0, "pt"),
  # zscore_legend_margin = margin(0, 0, 0, 10, "pt"),
  # ave_plot_margin = unit(c(80, 10, 0, 0), "pt"),
  # zscore_plot_margin = unit(c(80, -10, 0, 0), "pt"),
  # left_legend_pos = "left",
  # right_legend_pos = "right",
  # axis_title_x = element_blank(),
  # axis_title_y = element_blank(),
  # axis_title_x_lab = "",
  # axis_title_y_lab = "",
  # zscore_colors = viridis_pal(option = "cividis")(9), # cividis had a grey color in the middle that will be our zero
  # average_colors = viridis_pal(option = "viridis")(8), # viridis has no 0 grey/neutral color in the middle, which is good because this will be from zero to high
  # average_scaled_values = FALSE,
  # axis_text_y_remove_strings = c("KEGG_", "BIOCARTA_", "REACTOME_"),
  # axis_text_y_nchar = 25,
  # pdf_bookmark_title = "03_Clustered_Heatmap",
  # annotate_pdf = TRUE,
  # output_subfolder = "XX_Heatmap",
  # tile_width_factor = 1.01,
  # tile_height_factor = 1.01,
  
  
  t.colors=NA,
  fileName="cluster.cdt",
  linkage="complete",
  distance="pearson",
  contrast=2,
  returnSampleClust=F,
  rowNames=NA,
  rightMar=7,
  bottomMar=1
){
  library(viridis)
  my_script = "view_clustered_myHeatmap.R"
  my_annotation = paste0("Produce clustered heatmap view: ", my_script)
  a = function(new_text){
    env = parent.env(environment())
    assign('my_annotation', c(get('my_annotation', envir = env), new_text), envir = env)
  }
  
  if(my_pack %>% is_not_null){
    if(sample_names %>% is.null){sample_names = my_pack$sample_names}
    if(groups %>% is.null){groups = my_pack$groups}
    if(input_df %>% is.null){input_df = my_pack$input_df}
    if(base_file_name %>% is.null){base_file_name = paste0(my_pack$base_file_name, get_heatmap_ending())}
    if(base_title %>% is.null()){base_title = my_pack$base_title}
    if(output_dir %>% is.null){output_dir = my_pack$output_dir}
    if(pdf_metadata %>% is.null){pdf_metadata = my_pack$pdf_metadata}
    if(models_argument_list %>% is.null){models_argument_list = my_pack$argument_list}
    if(my_model %>% is.null){my_model = my_pack$my_model}
    if(imported_annotation %>% is.null){imported_annotation = my_pack$annotation}
  }
  
  if(is.null(pdf_bookmark_title) || is.na(pdf_bookmark_title)){
    pdf_bookmark_title = gsub("_", " ", base_file_name)
  }
  
  # if(is.null(pdf_bookmark_id) || is.na(pdf_bookmark_id)){
  #   pdf_bookmark_id = base_file_name
  # }

  a(paste0("Received data from ", my_model))
  
  a(format_imported_annotation(imported_annotation))
  a("")
  
  output_dir = file.path(output_dir, output_subfolder)
  dir.create(output_dir, showWarnings = FALSE)  # make output folder
  
  
  # we took care of all of the scaling from the unsupervised clustering 
  gene_matrix = t(as.matrix(input_df))
  colnames(gene_matrix) = sample_names



  library(ctc)
  
  cols <- function(lowi = "green", highi = "red", ncolors = 20) {
    low <- col2rgb(lowi)/255
    high <- col2rgb("black")/255
    col1 <- rgb(seq(low[1], high[1], len = ncolors), seq(low[2], 
                                                         high[2], len = ncolors), seq(low[3], high[3], len = ncolors))
    low <- col2rgb("black")/255
    high <- col2rgb(highi)/255
    col2 <- rgb(seq(low[1], high[1], len = ncolors), seq(low[2], 
                                                         high[2], len = ncolors), seq(low[3], high[3], len = ncolors))
    col<-c(col1[1:(ncolors-1)],col2)
    return(col)
  }
  
  
  temp<-hclust2treeview(
    gene_matrix,
    method=distance,
    file=file.path(output_dir, "cluster.cdt"),
    link=linkage,
    keep.hclust=T
    )
  gTree<-temp[[1]]
  sTree<-temp[[2]]
    
  imageVals<-gene_matrix
  imageVals[gene_matrix > contrast] <- contrast
  imageVals[gene_matrix < -1 * contrast] <- -1 * contrast
  
  if(sum(is.na(t.colors))>0){
    heatmap(imageVals,Rowv=as.dendrogram(gTree),Colv=as.dendrogram(sTree),
            col=cols(),labCol=colNames, scale="none",
            margins=c(bottomMar,rightMar),labRow=rowNames)
  }else{
    if(length(t.colors)>dim(imageVals)[2]){
      heatmap.plus(imageVals,Rowv=as.dendrogram(gTree),Colv=as.dendrogram(sTree),
                   col=cols(),labCol=colNames,labRow=rowNames,scale="none",
                   ColSideColors=t.colors, margins=c(bottomMar,rightMar))
    }else{
      heatmap(imageVals,Rowv=as.dendrogram(gTree),Colv=as.dendrogram(sTree),
              col=cols(),labCol=colNames,labRow=rowNames,scale="none",
              ColSideColors=as.vector(t(t.colors)), margins=c(bottomMar,rightMar))
    }
  }
  if(returnSampleClust){
    return(sTree)
  }
  
  write.table(output_df, file.path(output_dir, paste0(base_file_name, get_matrix_extension())), sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  heatmap_save_path = file.path(output_dir, paste0(base_file_name, ".pdf"))
  ggsave(heatmap_save_path, heat_maps, device = "pdf", width = pdf_width, height = pdf_height)
  
  if(annotate_pdf){
    # pdf_metadata is still being passed and could be added if we want, but for now just adding what we need for latex to make the bookmarks
    # format is a string c("value_id1: value1", "value_id=2: value2" 
    # to combine into nice pdfs we need the pdf_width, pdf_height, significance,	pdf_bookmark_title	pdf_bookmark_id

    # for significance if an fdr correction is shown we will output the highest number of stars
    
    my_meta_data = c(
      paste0("pdf_bookmark_title: ", pdf_bookmark_title),
      paste0("significance: ", max_significance)
    )
    
    write_metadata(heatmap_save_path, my_meta_data)
    
  }
  
  
  #write_metadata(heatmap_save_path, pdf_metadata)
  if(log_completion){
    warning("Haven't fixed toolkit's add_to_completed_analyses yet")
#     
#     add_to_completed_analyses(
#     classifier = models_argument_list$classifier,
#     included = models_argument_list$inclusion_list,
#     excluded = models_argument_list$exclusion_list,
#     model = my_model,
#     view = "Heatmap",
#     path = heatmap_save_path)
  }
  a('Heatmap Done' %>% as.footer)
  annotation_path = file.path(output_dir, paste0(base_file_name, get_note_extension()))
  write.table(my_annotation, annotation_path, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  
  return(list(my_path = heatmap_save_path, my_plot = heat_maps))
}
