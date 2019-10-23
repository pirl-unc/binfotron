# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# view_boxplot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#
# 
#
#' @title Create a boxplot of analysis results.
#'
#' @description
#' \code{view_boxplot} was created to make boxplots of any numeric features. pValues are calculated by AOV
#'   and displayed as stars over the significant features.  Red stars are fdr significant.
#'
#' These list items may be overridden by specifying them individually.  They are described in 
#' more detail below:
#' @param add_sample_counts Boolean to indicate whether the groupings should show the sample count next to them.
#' @param axis_text_y String sets the value for the y axis label.
#' @param boxplot_color String value goes into \code{geom_boxplot(color = boxplot_color)}. Good for overriding the scale_color_manual_values for the boxplot.
#' @param boxplot_fill String value goes into \code{geom_boxplot(fill = boxplot_fill)}.  Good for overriding the scale_fill_manual_values for the boxplot.
#' @param display_stars Boolean to indicate whether stars should be put over significant values
#' @param fdr_method Specifies the method for FDR correction. Typically use "none" or "BH". See \code{\link{p.adjust}} for other options.
#' @param feature_strip_position String indicating where the feature names should go.  'top', 'bottom', 'left' or 'right' Anything other than 'bottom' will likely need bebugging. \code{facet_wrap(~ feature, strip.position = feature_strip_position)}
#' @param feature_text_angle Adjusts angle of x-axis text: \code{strip.text.x = element_text(angle = feature_text_hjust)}
#' @param feature_text_btm_margin Adjusts bottom margin of x-axis text: \code{margin = margin(feature_text_top_margin, 0, feature_text_btm_margin, 0, "cm")}
#' @param feature_text_hjust Adjusts hjust of x-axis text: \code{strip.text.x = element_text(hjust = feature_text_hjust)}
#' @param feature_text_size Sets the size of the font for the x axis feature labels: \code{strip.text.x = element_text(size = feature_text_size)}
#' @param feature_text_top_margin Adjusts top margin of x-axis text: \code{margin = margin(feature_text_top_margin, 0, feature_text_btm_margin, 0, "cm")}
#' @param feature_text_vjust Adjusts vjust of x-axis text: \code{strip.text.x = element_text(vjust = feature_text_vjust)}
#' @param ggplot_theme theme element to make major changes to the plot appearance
#' @param graph_max Set max y-value of the graph.
#' @param graph_min Set min y-value of the graph.
#' @param legend_name String sets the name that will go next to the legend. \code{NULL} will use the grouping_col value. For no label, use double quotes ("").
#' @param my_position_dodge Sets the position_dodge for spacing the groups
#' @param my_subtitle String to specify the subtitle of the graph.
#' @param number_of_rows Sets how many rows of graphs will be made: \code{facet_wrap(~ feature, nrow = number_of_rows)}
#' @param order_features "ascending" "descending" or NULL
#' @param point_alpha Numeric to set the alpha value of the points on the graph.
#' @param point_color String value goes into \code{geom_point(color = point_color)}. Good for overriding the scale_fill_manual_values for the points. 
#' @param point_fill String value goes into \code{geom_point(fill = point_fill)}.
#' @param point_shape Sets the shape (pch) of the points. See http://sape.inf.usi.ch/quick-reference/ggplot2/shape for options.16-20 are good options for color only.  21-25 are point options with fills.
#' @param point_size Numeric to set the size of the points on the graph.
#' @param pValue_decimal_places Sets how many decimal places are displayed for pValues in the output stats file.
#' @param remove_feature_underscores Boolean to indicate whether the underscore should be removed from the features names.
#' @param scale_color_manual_values vector fo colors to use in coloring the points. Some examples:viridis_pal(alpha = 1, begin = 0, end = 1, direction = 1)(5); rainbow(10); or c("#FF0000FF" "#FF9900FF" "#CCFF00FF")
#' @param scale_fill_manual_values scale_fill_manual values. Some examples:viridis_pal(alpha = 1, begin = 0, end = 1, direction = 1)(5); rainbow(10); or c("#FF0000FF" "#FF9900FF" "#CCFF00FF")If there are fewer fill values then there are groups the values will be repeated to make enough. Important to note that not all point shapes have a 'fill'. See 'point_shape' for more details.
#' @param show_box_and_violin_legend boolean to indicate if box and violin legend should be shown
#' @param show_point_legend boolean to indicate if point legend should be shown
#' @param show_violins Boolean to indicate whether violins should be shown on the plot or not.  
#' @param star_text_size Set the font size of the significance stars
#' @param star_text_vjust Set the vjust of the significance stars
#' @param star_y_position Set the height of the bar for features with significant differences between the groups.
#' @param x_axis_line_size x_axis_line_size
#' @param show_points Boolean to indicate whether the the points should be shown.
#' @param outlier_color Color fills in outlier.color of geom_boxplot
#' @param outlier_size Integer for in outlier.size of geom_boxplot
#' @param violin_scale See \code{\link{geom_violin}} scale. "area" or "count"
#' @param violin_alpha Sets the alpha for the violins
#' @param violin_color String value goes into \code{geom_violin(color = violint_color)}. Good for overriding the scale_color_manual_values for the violin.
#' @param violin_line_width width of violin line
#' 
#' @return Several outputs are produced: 
#'   1) ggplot
#'   2) attribute(ggplot)$comments tells steps done by script in making the plot
#'   3) attribute(ggplot)$stats has a data.table with the stats calculated by the plot
#'   
#'     
#' @section Limitations:
#' \itemize{
#'   \item If samples are matched and one of them has an \code{NA} it won't drop it's match. 
#'   \item If samples are matches and there are different numbers of matches across conditions, it won't know what to do with this.
#' }
#' 
#' @section Future changes:
#' \itemize{
#'   \item Add option for doing aov or mann-whitney or Kolmogorov–Smirnov test
#'   \item Add option to drop samples with any missing features
#' }
#' 
#' @family view gene_signature
#' 
#' @export
view_boxplot = function(# passed from analysis app
  # Primary arguments
  base_title = NULL,
  feature_names = NULL,
  matched_col = NULL,
  input_dt = NULL,
  grouping_col = NULL,
  # ---
  add_sample_counts = TRUE,
  axis_text_y = "Gene Signature Expression",
  boxplot_color = NULL,
  boxplot_fill = NULL,
  boxplot_width = 0.9,
  display_stars = TRUE,
  fdr_method = "BH",
  feature_strip_position = "bottom",
  graph_max =  4,
  graph_min = -4,
  legend_name = NULL,
  my_position_dodge = 0.9,
  my_subtitle = NULL,
  number_of_rows = 2,
  order_features = NULL,
  point_alpha = 1,
  point_color = NULL,
  point_fill = NULL,
  point_shape = 20, 
  point_size = 0.2,
  pValue_decimal_places = 4,
  remove_feature_underscores = TRUE,
  scale_color_manual_values = "black",
  scale_fill_manual_values = gg_color_hue(5), 
  show_box_and_violin_legend = TRUE,
  show_point_legend = FALSE,
  show_violins = TRUE,
  star_text_size = 4,
  star_text_vjust = 0.4,
  star_y_position = graph_max*0.90,
  show_points = TRUE,
  outlier_color = NA,
  outlier_size = 1,
  violin_scale = "width",
  violin_alpha = 1,
  violin_color = NA,
  violin_line_width = 1,
  x_axis_line_size = 0,
  ggplot_theme = theme(
    legend.key.height=unit(1,"line"),
    legend.key.width=unit(1,"line"),
    legend.position = "top",
    legend.justification = "center",
    legend.key.size =  unit(0.1, "in"),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 12),
    plot.subtitle = element_text(size=8, hjust=0.5, face="italic"),
    plot.title = element_text(size=18, face='bold', hjust= 0.5),
    legend.text = element_text(size=9),
    legend.title = element_text(size=11),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_line(size = 0.5),
    axis.ticks.y = element_line(size = 0.5),
    panel.spacing.x = unit(2, "mm"),
    panel.spacing.y = unit(0, "mm"),
    strip.background = element_rect(fill = NA),
    strip.text.x = element_text(
      size = 8, 
      hjust = 1, 
      vjust = 0.5,
      margin = margin(0.1, 0, 0.4, 0, "cm"), 
      angle = 90
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "lightgrey", size=0.5),
    panel.background = element_rect(fill = "white",
                                    colour = "white")
  )
){
  
  library(plyr)
  library(tidyr)
  
  my_script = "view_boxplot.R"
  my_annotation = paste0("Produced boxplot view: ", my_script)
  a = function(new_text){
    env = parent.env(environment())
    assign('my_annotation', c(get('my_annotation', envir = env), new_text), envir = env)
  }
  
  if("comments" %in% names(attributes(input_dt))){
    a("Input data comments:")
    a(format_imported_annotation(attributes(input_dt)$comments))
  }

  boxplot_df = input_dt %>% as.data.frame()
  boxplot_df = boxplot_df[complete.cases(boxplot_df[[grouping_col]]),] # to get group counts
  boxplot_df[[grouping_col]] = factor(boxplot_df[[grouping_col]])
  # get the n's 
  category_n = summary(boxplot_df[[grouping_col]])


  boxplot_df[,feature_names ] = apply(boxplot_df[,feature_names, drop = FALSE], 2, function(x){as.numeric(as.character(x))}) # need to make the columns the same type or we will get a complaint that data is lost in 
  if(!is.null(matched_col)){
    
    boxplot_df = boxplot_df[, c(matched_col, grouping_col, feature_names)] # this converts pipes into periods in the titles so names are fixed in next step
    names(boxplot_df) = c("Patient_ID", "sample_groups", names(input_df))
    
    boxplot_df = gather(boxplot_df, key = grouping_col, value=feature_names, 3:ncol(boxplot_df))
    names(boxplot_df) = c("Patient_ID", grouping_col, "feature", 'value')
    
  } else {
    
    boxplot_df = boxplot_df[, c(grouping_col, feature_names)] # this converts pipes into periods in the titles so names are fixed in next step
    boxplot_df = boxplot_df[complete.cases(boxplot_df),]
    
    boxplot_df = gather(boxplot_df, key = grouping_col, value=feature_names, 2:ncol(boxplot_df))
    names(boxplot_df) = c(grouping_col, "feature", 'value')
    
  }
  
  boxplot_df = boxplot_df[complete.cases(boxplot_df),]
  
  boxplot_df$feature = factor(boxplot_df$feature, levels = feature_names)
  
  if(!is.null(order_features)){
    # first get 50th percentile for each feature
    fiftieth_percentile = c()
    new_order = NULL
    for (my_feature in levels(boxplot_df$feature)){
      subdat = boxplot_df[boxplot_df$feature == my_feature, ]
      fiftieth_percentile = c(fiftieth_percentile, median(subdat$value, na.rm = TRUE))
    }
    if (order_features == "ascending"){
      
      new_order = order(fiftieth_percentile, decreasing = FALSE)
      
    } else if((order_features == "descending")) {
      
      new_order = order(fiftieth_percentile, decreasing = TRUE)
      
    } 
    if (!is.null(new_order)) {
      boxplot_df$feature = factor(boxplot_df$feature, levels(boxplot_df$feature)[new_order])
    }
  }

  # Do stats on the data ----------------------------------------------------------------------
  output_feature = vector()
  output_aov = vector()
  f_value = vector()
  fdr_aov = vector()
  output_n = vector()
  
  if(matched_col %>% is.null){
    a("Calculated error using AOV.  Samples were not matched.")
  } else {
    a("Calculated error using AOV accounting for matched samples.")
  }
  a("")
  
  for (feature in feature_names) {
    temp_df = boxplot_df[boxplot_df$feature == feature, ]
    temp_df = temp_df[complete.cases(temp_df$value), ]
    output_feature = c(output_feature, feature)
    
    if(length(levels(factor(temp_df[[grouping_col]]))) < 2){
      output_aov = c(output_aov, NA)
      f_value = c(f_value, NA)
      feature_msg = paste0("The feature ", feature, " was not present in more than 1 condition so stats could not be done.")
      warning(feature_msg)
      a(feature_msg %>% wrap_sentence)
    } else {
      
      if(!is.null(matched_col)){
        # some matches might have an NA for this feature.  need to get rid of those
        item_counts = summary(temp_df[["Patient_ID"]]) 
        lone_items = names(item_counts[item_counts < 2])
        if(length(lone_items) > 0){
          temp_df = temp_df[temp_df[["Patient_ID"]] %ni% lone_items, ]
        }
        
        my_test = summary(aov(temp_df$value~temp_df[[grouping_col]] + Error(temp_df[["Patient_ID"]])))
        output_aov = c(output_aov, unlist(my_test[[2]])[['Pr(>F)1']])
        f_value = c(f_value, NA)
        
      } else {
        f_value = c(f_value, summary(aov(temp_df$value~temp_df[[grouping_col]]))[[1]][["F value"]][1])
        output_aov = c(output_aov, summary(aov(temp_df$value~temp_df[[grouping_col]]))[[1]][["Pr(>F)"]][1])
      }
    }
  }

  
  
  
  
  # FDR correct -----------------------------------------------------------------------------------
  max_significance = NA
  if (!is.null(fdr_method)) {
    a(paste0("FDR corrected using ", fdr_method, "."))
    fdr_aov = p.adjust(output_aov, method= fdr_method) # Benjamini & Hochberg aka. FDR
    pValue_label = paste0(fdr_method,'_pValue')
    if(!all(is.na(fdr_aov))){
      max_significance = min(fdr_aov, na.rm = T)
    }
  } else {
    a(paste0("No FDR correction was made."))
    pValue_label = "pValue"
    if(!all(is.na(output_aov))){
      max_significance = min(output_aov, na.rm = T)
    }
  }
  
  if(!is.na(max_significance)) {
    max_significance = pvalue_stars(max_significance)
  } else {
    warning("The max_significance was NA. This could mean your category only had one level in it. Does it?")
    display_stars = FALSE
  }
  
  a("")
  
  if(add_sample_counts){
    for(my_level in levels(boxplot_df[,grouping_col])){
      my_count = category_n[my_level] %>% as.numeric()
      boxplot_df[,grouping_col] = mapvalues(boxplot_df[,grouping_col], from = c(my_level), to = paste0(my_level, "\n(n=",my_count,")"))
    }
  }
  
  
  output_stats_df = data.frame(output_feature, output_aov, fdr_aov, f_value)
  names(output_stats_df) = c('feature', 'pValue', pValue_label, "F_Value")
  
  
  if(remove_feature_underscores){
    for(my_level in levels(boxplot_df[,"feature"])){
      boxplot_df[,"feature"] = mapvalues(boxplot_df[,"feature"], from = c(my_level), to = gsub("_", " ", my_level))
    }
    for(my_level in levels(output_stats_df[,"feature"])){
      output_stats_df[,"feature"] = mapvalues(output_stats_df[,"feature"], from = c(my_level), to = gsub("_", " ", my_level))
    }
  }
  
  my_groups = levels(boxplot_df[ , grouping_col])
  
  enough_scale_color_manual_values = scale_color_manual_values
  if(length(enough_scale_color_manual_values) < length(my_groups)){
    enough_scale_color_manual_values = rep.int(enough_scale_color_manual_values, length(my_groups))[1:length(my_groups)]
  }
  
  enough_scale_fill_manual_values = scale_fill_manual_values
  if(length(enough_scale_fill_manual_values) < length(my_groups)){
    enough_scale_fill_manual_values = rep.int(enough_scale_fill_manual_values, length(my_groups))[1:length(my_groups)]
  }
  
  
  geom_point_args = c(
    "alpha = point_alpha", 
    "position = position_dodge(my_position_dodge)", 
    "size = point_size", 
    "show.legend = show_point_legend", 
    "shape = point_shape"
  )
  
  if(!is.null(point_color)){
    geom_point_args = c(geom_point_args, "color = point_color")
  }
  
  if(!is.null(point_fill)){
    geom_point_args = c(geom_point_args, "fill = point_fill")
  }
  
  
  geom_boxplot_args = c(
    "aes_string(x = grouping_col, y = 'value', fill = grouping_col)",
    "size = 0.3", # box line width
    paste0("outlier.colour = '", outlier_color,"'"), 
    paste0("outlier.size = ", outlier_size),
    "show.legend = show_box_and_violin_legend",
    "position = position_dodge(my_position_dodge)", 
    "width = boxplot_width"
  )
  
  if(!is.null(boxplot_color)){
    geom_boxplot_args = c(geom_boxplot_args, "color = boxplot_color")
  }
  
  if(!is.null(boxplot_fill)){
      geom_boxplot_args = c(geom_boxplot_args, "fill = boxplot_fill")
  } else {
    if( show_violins ){
      geom_boxplot_args = c(geom_boxplot_args, "fill = NA") # it's a little tacky to show the boxplot fills with the violins
    } 
  }
  
  
  my_plot = ggplot(boxplot_df, aes_string(x = grouping_col, y = "value", fill = grouping_col,  color = grouping_col))
  if( show_violins ){
    if(is.null(violin_color)){
      my_plot = my_plot + 
        geom_violin(alpha = violin_alpha, show.legend = show_box_and_violin_legend, scale = violin_scale, size = violin_line_width)
    } else {
      my_plot = my_plot + 
        geom_violin(alpha = violin_alpha, color = violin_color, show.legend = show_box_and_violin_legend, scale = violin_scale, size = violin_line_width)
      
    }
  } 
  my_plot = my_plot +
    eval(parse(text = paste0("geom_boxplot(", paste0(geom_boxplot_args, collapse = ", "), ")")))
  
  if(show_points){
    my_plot = my_plot +
    eval(parse(text = paste0("geom_point(", paste0(geom_point_args, collapse = ", "), ")")))
  }
  
  if(display_stars){
    # have to set this up before the plot so we'll know if the subtitle needs to talk about 
    #  red or black stars
    
    library(ggsignif)
    
    show_stars_df = output_stats_df[output_stats_df$pValue <= 0.05 | output_stats_df[,pValue_label] <= 0.05, ]

    if(nrow(show_stars_df) > 0){
      show_stars_df$Stars = NA
      show_stars_df$Color = NA
      my_groups = levels(boxplot_df[,grouping_col])
      show_stars_df$Start = my_groups[1]
      show_stars_df$End = my_groups[length(my_groups)]
      for(row_index in 1:nrow(show_stars_df)){
        if(show_stars_df[row_index, pValue_label] <= 0.05){
          show_stars_df$Stars[row_index] = pvalue_stars(show_stars_df[row_index, pValue_label])
          show_stars_df$Color[row_index] = "red"
        } else if (show_stars_df[row_index, "pValue"] <= 0.05) {
          show_stars_df$Stars[row_index] = pvalue_stars(show_stars_df[row_index, "pValue"])
          show_stars_df$Color[row_index] = "black"
        }
      }
      
    } else {
      display_stars = FALSE
    }
    
    if(is.null(my_subtitle)){
      
      if(any(show_stars_df$Color == "red")){
        my_subtitle = "(red stars - FDR significant)"
      } else if(any(show_stars_df$Color == "black")){
        my_subtitle = "(black stars - not FDR significant)"
      } else {
        my_subtitle = "(none significant)"
      }
    }
  }
  
  my_plot = my_plot +
    scale_color_manual(values = enough_scale_color_manual_values) + 
    scale_fill_manual(values = enough_scale_fill_manual_values) +
    facet_wrap(~ feature, nrow = number_of_rows, strip.position = feature_strip_position) + 
    labs(title = base_title,  x = NULL, y = axis_text_y,
         subtitle = my_subtitle, fill = legend_name, color = legend_name) + 
    coord_cartesian(ylim = c(graph_min, graph_max), expand = F) +
    geom_hline(yintercept = graph_min, size = x_axis_line_size) +
    ggplot_theme


  # color here is buggy.  they want 3 colors for each label: 
  #   the first is the label and the left tip, 
  #   the second is the long bar, 
  #   the third is the right tip
  # since we don't want one tip to be the same color as the label, we set the tip length to 0.00
  #   (setting to NA makes it count the colors differently)
      
  if(display_stars){
    stat_colors =  unlist(lapply(show_stars_df$Color, function(x){c(x,"black", "black")})) # sapply would have returned a matrix here
    my_plot = my_plot +
      ggsignif::geom_signif(data=show_stars_df, inherit.aes = FALSE,
                  mapping = aes_string(xmin="Start", xmax="End", annotations="Stars"), y_position=star_y_position, color = stat_colors, # the three indicates different parts of the sig marker
                  textsize = star_text_size, vjust = star_text_vjust, tip_length = 0.00,
                  manual=TRUE)
    cat("geom_signif: gives a fake warning about 'Ignoring unknown aesthetics'. This warning should be ingnored.\nFor more info see https://cran.r-project.org/web/packages/ggsignif/README.html\n")
  }
  
  output_stats_df[,"pValue"] = specify_decimal(output_stats_df[,"pValue"], pValue_decimal_places)
  output_stats_df[,pValue_label] = specify_decimal(output_stats_df[,pValue_label], pValue_decimal_places)
  
  attributes(my_plot)$stats = output_stats_df
  
  attributes(my_plot)$comments = my_annotation
  
  return(my_plot)
  
}
