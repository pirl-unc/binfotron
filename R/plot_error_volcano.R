# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot_error_volcano
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Creates a scatterplot with optional errorbars and fine control over labels, colors, alpha and more
#' 
#' @param plot_df Dataframe containing at a minimum x and y axis data
#' @param axis_clm_x column name in plot_df corresponding to x axis
#' @param axis_clm_y column name in plot_df corresponding to y axis
#' 
#' @param alpha_clm column name of categorical or continuous variable in plot_df on which to base alpha aesthetic
#' @param color_clm column name of categorical variable in plot_df on which to base color aesthetic
#' @param color_scheme ggplot function to specify discrete colors
#' @param error_clm_lower column name containing lower bound of error bars
#' @param error_clm_upper column name containing upper bound of error bars
#' @param gg_add_on ggplot elements to add on to plot. Defaults to \code{theme_minimal()}
#' @param independent_clm column name in plot_df from which to pull values for reporting on removed extreme values. If not set, and extreme values are removed, only the number of removed observations will be reported.
#' @param lab_color String for labeling color legend
#' @param lab_x String for labeling x axis
#' @param lab_y String for labeling y axis
#' @param label_clm column name in plot_df from which to pull labels for points
#' @param label_size numeric value to use for the size of the point labels
#' @param max_labels Integer to set the max number of labels that will be shown
#' @param output_path Output path for plot
#' @param plot_title String for plot title
#' @param size_height Numeric to specify plot height
#' @param size_width Numeric to specify plot width
#'  
#' @return Returns ggplot object with point and optionally errorbar geoms as well as point labels if desired
#' 
#' @export
#' 
plot_error_volcano = function( 
	plot_df, 
	axis_clm_x,
	axis_clm_y,
	alpha_clm=NA,
	color_clm=NA,
	color_scheme=scale_color_viridis_d(option="rocket",begin = 0.8, end=0.1),
	error_clm_lower = NA,
	error_clm_upper = NA,
	gg_add_on = theme_minimal(),
	independent_clm=NA,
	lab_color="Q-Value",
	lab_x="Hazard Ratio (Log10)",
	lab_y="Uncorrected Significance (-Log10)",
	label_clm=NA,
	label_nudge_x = 0.025,
	label_nudge_y = -0.15,
	label_size=1.5,
	max_labels = 20,
	output_path=NA,
	plot_title="Significance vs. Hazard Ratio",
	size_height = 4.5,
	size_width = 6
){
	#for easier reference down the line, get colnames from plot_df
	dtp_clms <- colnames(plot_df)
	
	# clear labels where alpha == 0 
	if(alpha_clm %in% dtp_clms & label_clm %in% dtp_clms){
		plot_df[[label_clm]][plot_df[[alpha_clm]]<1] = ""
	}
	
	#ggplot throws an error if NA is passed to the color aesthetic ( although not to alpha or label ... ), use NULL for default instead
	
	#create base plot with axes and scatterplotted data
	base_plot = ggplot(plot_df,aes_string(axis_clm_x,axis_clm_y, 
																				label=label_clm)) +
		geom_vline(xintercept=0, size=0.25)+
		geom_hline(yintercept=0, size=0.25)+
		geom_point(mapping = aes_string(color=color_clm), alpha = plot_df[[alpha_clm]])
	
	#if errorbars were requested, add them to the plot
	if( !is.na(error_clm_lower)  && !is.na(error_clm_upper) ){
		#two values must be defined as upper and lower error columns
		if (error_clm_lower %in% dtp_clms & error_clm_upper %in% dtp_clms){
			base_plot <- base_plot + 
				geom_errorbarh(
					mapping = aes_string(xmin=error_clm_lower,xmax=error_clm_upper, color=color_clm),
					alpha = plot_df[[alpha_clm]]
				)
		}
	}
	
	label_df = plot_df
	if (sum(plot_df[[label_clm]] != "") > max_labels){
		label_df = plot_df[order(plot_df[[axis_clm_y]], decreasing = T),]
		label_df = label_df[1:max_labels, ]
	}
	
	#if labels were included in the plot, use ggrepel to format them properly
	if( is.character(label_clm) ){
		base_plot <- base_plot + 
			geom_label_repel(
				data = label_df,
				show.legend = F,
				size = label_size,
				label.padding = 0.1,
				force = 1,
				direction='both',
				nudge_x = label_nudge_x,
				nudge_y = label_nudge_y,
				min.segment.length = 0,
				segment.size = 0.25
			)
	}
	
	base_plot <- base_plot + labs(title=plot_title,x=lab_x,y=lab_y, color=lab_color)+
		theme(plot.title = element_text(hjust = 0.5))
	
	#if either color, alpha or both are defined for this plot, hide the alpha legend and use the inferno color scale to make things pretty
	if( is.character(color_clm) || is.character(alpha_clm) ){
		base_plot <- base_plot + guides(alpha="none") +
			color_scheme
	}
	
	if ( !is.null(gg_add_on) ) base_plot = base_plot + gg_add_on
	
	#if the graph is to be printed to file, do so now
	if(!missing(output_path)) ggsave(output_path, base_plot, width = size_width, height = size_height)
	
	if (rstudioapi::isAvailable()) base_plot
}
