# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot_cumulative_features
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Plots a line depicting the number of features ( or any categorical variable ) that are below the pvalue ( or any categorical variable )
#' 
#' @description Intended to depict how many features are below pvalues with vertical lines representing the typical cutoffs for significance ( .05, .01, .001 )
#' 
#' 
#' @param plot_df Dataframe containing at a minimum pvalue data
#' @param lab_x String for labeling x axis
#' @param lab_y String for labeling y axis
#' @param max_features integer representing the maximum number of features to plot. Values above this will be dropped from the plot_df
#' @param max_pvalue double representing the highest pvalue to be plotted. Values above this will be dropped from the plot_df
#' @param output_path Output path for plot
#' @param plot_title String for plot title
#' @param pvalue_clm column name in plot_df corresponding to pvalues
#' @param size_height Numeric to specify plot height
#' @param size_width Numeric to specify plot width
#' 
#' @return Returns ggplot object with graphed line and vertical lines at significance intervals
#' 
#' @export
#' 
plot_cumulative_features <- function( 
	plot_df, 
	lab_x="P-Value ( Log10 )",
	lab_y="Cumulative Feature Count ( Log10 )",
	max_features=NA,
	max_pvalue=NA,
	output_path=NA,
	plot_title="Cumulative Features Below P-Value",
	pvalue_clm="pValue",
	size_height = 5,
	size_width = 5.5
){
	
	# does pvalue_clm exist in data?
	if( !(pvalue_clm %in% names(plot_df)) ) stop("Value sent for pvalue_clm (", pvalue_clm, ") does not exist in data.")
	
	# order by ascending pvalue column
	plot_df %<>% .[order(.[[pvalue_clm]]),]
	
	# add a 'Features_Below' field which represents the # of features below each pValue
	plot_df$Features_Below <- c(1:nrow(plot_df))
	
	# if a max_pvalue was sent, remove all observations below that threshold
	if( !is.na(max_pvalue) ){
		plot_df %<>% .[ .[ pvalue_clm ] < max_pvalue, ]
	}
	
	# if a max_features was sent, remove all observations below that threshold
	if( !is.na(max_features) ){
		plot_df %<>% .[ .$Features_Below <= max_features, ]
	}
	
	#plot it at last
	base_plot <- ggplot(plot_df, aes_string(x=pvalue_clm, y="Features_Below")) + # pValue vs. count below pValue
		geom_line() + #color="blue") + # yep, plot a line ...
		scale_x_log10(labels=function(x) format(x, scientific = TRUE), breaks = scales::trans_breaks("log10", function(x) 10^x)) +# scale x axis on Log10
		scale_y_log10() + #also scale y axis on Log10
		geom_vline(xintercept=c(.001,.01,.05), color="indianred1", linetype="dashed") + # add some vertical lines at important significance levels
		annotation_logticks(side="bl") + # add tickmarks along the bottom for the log scale
		labs(title=plot_title, x=lab_x, y=lab_y) + # add labels for title, x and y axis
		theme(plot.title = element_text(hjust = 0.5), # center title
					panel.grid.minor.x = element_blank(), # remove minor x gridlines
					panel.grid.major.x = element_blank()) + # remove major x gridlines
		annotate("text", x=.0008, y=nrow(plot_df), label="0.001", angle=90) + #label the significance level vertical lines ... this is pretty janky
		annotate("text", x=.008, y=nrow(plot_df), label="0.01", angle=90) +
		annotate("text", x=.0405, y=nrow(plot_df), label="0.05", angle=90) #+ 
	#annotate("text", x=.16, y=nrow(plot_df), label="0.20", angle=90)
	
	#print it to file if desired
	if(!missing(output_path)) ggsave(output_path, base_plot, width = size_width, height = size_height)
	
	if (rstudioapi::isAvailable()) base_plot
}

