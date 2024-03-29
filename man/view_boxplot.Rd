% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/view_boxplot.R
\name{view_boxplot}
\alias{view_boxplot}
\title{Create a boxplot of analysis results.}
\usage{
view_boxplot(
  base_title = NULL,
  feature_names = NULL,
  matched_col = NULL,
  input_dt = NULL,
  grouping_col = NULL,
  add_sample_counts = TRUE,
  axis_text_y = "Gene Signature Expression",
  boxplot_color = NULL,
  boxplot_fill = NULL,
  boxplot_width = 0.9,
  display_stars = TRUE,
  fdr_method = "BH",
  feature_strip_position = "bottom",
  graph_max = 4,
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
  star_y_position = graph_max * 0.9,
  show_points = TRUE,
  outlier_color = NA,
  outlier_size = 1,
  violin_scale = "width",
  violin_alpha = 1,
  violin_color = NA,
  violin_line_width = 1,
  x_axis_line_size = 0,
  ggplot_theme = theme(legend.key.height = unit(1, "line"), legend.key.width = unit(1,
    "line"), legend.position = "top", legend.justification = "center", legend.key.size =
    unit(0.1, "in"), axis.text.y = element_text(size = 9), axis.title.y =
    element_text(size = 12), plot.subtitle = element_text(size = 8, hjust = 0.5, face =
    "italic"), plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 9), legend.title = element_text(size = 11),
    axis.text.x = element_blank(), axis.line.x = element_blank(), 
     axis.ticks.x =
    element_blank(), axis.line.y = element_line(size = 0.5), axis.ticks.y =
    element_line(size = 0.5), panel.spacing.x = unit(2, "mm"), panel.spacing.y = unit(0,
    "mm"), strip.background = element_rect(fill = NA), strip.text.x = element_text(size =
    8, hjust = 1, vjust = 0.5, margin = margin(0.1, 0, 0.4, 0, "cm"), angle = 90),
    panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "lightgrey", size = 0.5), panel.background
    = element_rect(fill = "white", 
         colour = "white"))
)
}
\arguments{
\item{add_sample_counts}{Boolean to indicate whether the groupings should show the sample count next to them.}

\item{axis_text_y}{String sets the value for the y axis label.}

\item{boxplot_color}{String value goes into \code{geom_boxplot(color = boxplot_color)}. Good for overriding the scale_color_manual_values for the boxplot.}

\item{boxplot_fill}{String value goes into \code{geom_boxplot(fill = boxplot_fill)}.  Good for overriding the scale_fill_manual_values for the boxplot.}

\item{display_stars}{Boolean to indicate whether stars should be put over significant values}

\item{fdr_method}{Specifies the method for FDR correction. Typically use "none" or "BH". See \code{\link{p.adjust}} for other options.}

\item{feature_strip_position}{String indicating where the feature names should go.  'top', 'bottom', 'left' or 'right' Anything other than 'bottom' will likely need bebugging. \code{facet_wrap(~ feature, strip.position = feature_strip_position)}}

\item{graph_max}{Set max y-value of the graph.}

\item{graph_min}{Set min y-value of the graph.}

\item{legend_name}{String sets the name that will go next to the legend. \code{NULL} will use the grouping_col value. For no label, use double quotes ("").}

\item{my_position_dodge}{Sets the position_dodge for spacing the groups}

\item{my_subtitle}{String to specify the subtitle of the graph.}

\item{number_of_rows}{Sets how many rows of graphs will be made: \code{facet_wrap(~ feature, nrow = number_of_rows)}}

\item{order_features}{"ascending" "descending" or NULL}

\item{point_alpha}{Numeric to set the alpha value of the points on the graph.}

\item{point_color}{String value goes into \code{geom_point(color = point_color)}. Good for overriding the scale_fill_manual_values for the points.}

\item{point_fill}{String value goes into \code{geom_point(fill = point_fill)}.}

\item{point_shape}{Sets the shape (pch) of the points. See http://sape.inf.usi.ch/quick-reference/ggplot2/shape for options.16-20 are good options for color only.  21-25 are point options with fills.}

\item{point_size}{Numeric to set the size of the points on the graph.}

\item{pValue_decimal_places}{Sets how many decimal places are displayed for pValues in the output stats file.}

\item{remove_feature_underscores}{Boolean to indicate whether the underscore should be removed from the features names.}

\item{scale_color_manual_values}{vector fo colors to use in coloring the points. Some examples:viridis_pal(alpha = 1, begin = 0, end = 1, direction = 1)(5); rainbow(10); or c("#FF0000FF" "#FF9900FF" "#CCFF00FF")}

\item{scale_fill_manual_values}{scale_fill_manual values. Some examples:viridis_pal(alpha = 1, begin = 0, end = 1, direction = 1)(5); rainbow(10); or c("#FF0000FF" "#FF9900FF" "#CCFF00FF")If there are fewer fill values then there are groups the values will be repeated to make enough. Important to note that not all point shapes have a 'fill'. See 'point_shape' for more details.}

\item{show_box_and_violin_legend}{boolean to indicate if box and violin legend should be shown}

\item{show_point_legend}{boolean to indicate if point legend should be shown}

\item{show_violins}{Boolean to indicate whether violins should be shown on the plot or not.}

\item{star_text_size}{Set the font size of the significance stars}

\item{star_text_vjust}{Set the vjust of the significance stars}

\item{star_y_position}{Set the height of the bar for features with significant differences between the groups.}

\item{show_points}{Boolean to indicate whether the the points should be shown.}

\item{outlier_color}{Color fills in outlier.color of geom_boxplot}

\item{outlier_size}{Integer for in outlier.size of geom_boxplot}

\item{violin_scale}{See \code{\link{geom_violin}} scale. "area" or "count"}

\item{violin_alpha}{Sets the alpha for the violins}

\item{violin_color}{String value goes into \code{geom_violin(color = violint_color)}. Good for overriding the scale_color_manual_values for the violin.}

\item{violin_line_width}{width of violin line}

\item{x_axis_line_size}{x_axis_line_size}

\item{ggplot_theme}{theme element to make major changes to the plot appearance}

\item{feature_text_angle}{Adjusts angle of x-axis text: \code{strip.text.x = element_text(angle = feature_text_hjust)}}

\item{feature_text_btm_margin}{Adjusts bottom margin of x-axis text: \code{margin = margin(feature_text_top_margin, 0, feature_text_btm_margin, 0, "cm")}}

\item{feature_text_hjust}{Adjusts hjust of x-axis text: \code{strip.text.x = element_text(hjust = feature_text_hjust)}}

\item{feature_text_size}{Sets the size of the font for the x axis feature labels: \code{strip.text.x = element_text(size = feature_text_size)}}

\item{feature_text_top_margin}{Adjusts top margin of x-axis text: \code{margin = margin(feature_text_top_margin, 0, feature_text_btm_margin, 0, "cm")}}

\item{feature_text_vjust}{Adjusts vjust of x-axis text: \code{strip.text.x = element_text(vjust = feature_text_vjust)}}
}
\value{
Several outputs are produced: 
  1) ggplot
  2) attribute(ggplot)$comments tells steps done by script in making the plot
  3) attribute(ggplot)$stats has a data.table with the stats calculated by the plot
}
\description{
\code{view_boxplot} was created to make boxplots of any numeric features. pValues are calculated by AOV
  and displayed as stars over the significant features.  Red stars are fdr significant.

These list items may be overridden by specifying them individually.  They are described in 
more detail below:
}
\section{Limitations}{

\itemize{
  \item If samples are matched and one of them has an \code{NA} it won't drop it's match. 
  \item If samples are matches and there are different numbers of matches across conditions, it won't know what to do with this.
}
}

\section{Future changes}{

\itemize{
  \item Add option for doing aov or mann-whitney or Kolmogorov–Smirnov test
  \item Add option to drop samples with any missing features
}
}

\concept{view gene_signature}
