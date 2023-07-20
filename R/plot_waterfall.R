# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot_waterfall
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Plots a bar graph depicting the hazard ratio ( on log10 scale ) by feature
#' 
#' @description Intended to aid in selection of significant features by displaying the hazard ratios for each feature ( ordered by hazard ratio ) with optional color coding for various Q-Value significance levels
#' 
#' @param plot_df Dataframe containing at a minimum hazard ratios and features
#' @param colors_clm column name in plot_df of categorical variable to use for coloring of bars ( typically representing significance levels )
#' @param features_clm column name in plot_df corresponding to features
#' @param fill_scheme ggplot function to specify discrete fills
#' @param hr_clm column name in plot_df corresponding to hazard ratios
#' @param lab_color String for labeling color legend
#' @param lab_x String for labeling x axis
#' @param lab_y String for labeling y axis
#' @param output_path Output path for plot
#' @param plot_title String for plot title
#' @param sig_clm column name in plot_df corresponding to q-values, only used if sig_threshold is also set
#' @param sig_threshold double representing the maximum q-value to plot on the graph. Observations with higher qvalues will be dropped from the data
#' @param size_height Numeric to specify plot height
#' @param size_width Numeric to specify plot width
#'  
#' @return Returns ggplot object with graphed line and vertical lines at significance intervals
#' 
#' @export
plot_waterfall <- function( 
	plot_df, 
	colors_clm=NA,
	features_clm="Independant",
	fill_scheme=scale_fill_viridis_d(option="rocket",begin = 0.8, end=0.3), 
	hr_clm="Hazard_Ratio",
	lab_color="Q-Value",
	lab_x="Features",
	lab_y="Hazard Ratio ( Log10 )",
	output_path=NA,
	plot_title=title_prefix,
	sig_clm="qValue",
	sig_threshold=NA,
	size_height = 5.5,
	size_width = 5,
	log10_scale = TRUE
){
	
	# do hr_clm and features_clm exist in data?
	if( !all(c(hr_clm, features_clm) %in% names(plot_df)) ) stop("Values sent for hr_clm and features_clm must exist in data.")
	
	# if a qval threshold is set, remove data values above that threshold
	if (sum(plot_df[[sig_clm]]<= sig_threshold) == 0) {
		sig_threshold = NA
		message(paste0("No features were below sig_threshold of ", sig_threshold, " so the threshold was dropped."))
	}
	
	if( is.numeric(sig_threshold) & sig_threshold > 0 & sig_threshold <= 1 ){
		if( !(sig_clm %in% names(plot_df))) stop("A sig_threshold was sent but the sig_clm does not exist in the data.")
		plot_df %<>% .[.[[sig_clm]] <= sig_threshold, ]
	}
	
	# subset data to be plotted to get around issue with ggplot aes evaluating parameters as character strings rather than column names
	df_to_plot <- data.frame(y=plot_df[[hr_clm]], x=plot_df[[features_clm]])
	colors_col = NULL
	
	# if colors are desired, ensure the proper data has been sent and set the colors_col data
	if(!is.na(colors_clm)){
		if(!(colors_clm %in% names(plot_df))) stop("Value sent for colors_clm must exist in data.")
		colors_col <- plot_df[[colors_clm]]
	}
	
	#plot it, finally
	if(log10_scale) {
		base_plot <- ggplot(df_to_plot, aes(reorder(x,y),log10(y)))
	} else {
		base_plot <- ggplot(df_to_plot, aes(reorder(x,y),y))
	}
	
	base_plot <- base_plot + # pValue vs. count below pValue
		geom_col(aes(fill=colors_col), width=1) + # yep, plot some bars ...
		geom_vline(xintercept=0) + # add y axis
		scale_x_discrete(breaks=NULL) + # get rid of x axis labels ... too many features to write them out as labels
		labs(title=plot_title, x=lab_x, y=lab_y, fill=lab_color) + # add labels for title, x and y axis
		fill_scheme + #apply nicer color scheme
		theme(plot.title = element_text(hjust = 0.5)) # center title
	
	#print it to file if desired
	if(!missing(output_path)) ggsave(output_path, base_plot, width = size_width, height = size_height)
	
	if (rstudioapi::isAvailable()) base_plot
}
