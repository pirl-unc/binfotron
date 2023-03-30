#library(data.table)
#library(magrittr)
#library(ggplot2)
#library(ggrepel)
#library(ggforce) # for geom_mark_ellipse

#' function dependencies:
#' 
#' find_central_elements_by_cluster
#'     get_ranks_from_df
#'     compare_via_mhorn
#'     isolate_central_cluster_elements
#'         pca_central_element
#'         correlation_central_element
#'     cluster_list_to_df
#'     write_list_as_gmt
#'     write_central_elements_table
#'     write_ranked_central_elements_table
#'     pca_variance_plot
#'     cluster_group_plot
#'     pc1_vs_pc2_plot
#'     
#' group_clusters_by_size - useful for looking at results
#' collapse_identical_clusters - useful for looking at results
#' print_cluster_counts: takes output from collapse_identical_clusters


get_ranks_from_df = function( ranked_df, rank_clm=NA ){
  if (is.na(rank_clm) | !(rank_clm %in% colnames(ranked_df))){
    cat("No rank_clm was provided, or rank_clm does not exist in data.frame. Proceeding with alphabetical ranking.")
    # sort alphabetically
    rtn <- c(1:nrow(ranked_df))
    names(rtn) <- sort(rownames(ranked_df))
  } else {
    rtn <- ranked_df[[rank_clm]]
    names(rtn) <- rownames(ranked_df)
    rtn %<>% sort()
  }
	
  return( rtn )
}


# Notes on find_central_elements_by_cluster:
# Error in kmeans(elements_df, cluster_ct, iter.max = x, nstart = y) : 
# 	more cluster centers than distinct data points.
# means data have fewer distinct cases than the number of centers specified. specify a lower max_clusters.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# find_central_elements_by_cluster
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Encapsulation of steps to create clusters and determine most central elements of each cluster
#' 
#' @description Generate clusters using kmeans method, and determine most representative element for each cluster using a pca analysis (most central feature in pca space) , mhorn similarity index (most similar feature), or pearson/spearman correlation (most correlated feature).
#' 
#' 
#' @param feature_df data.frame on which to perform PCA, mhorn or spearman analysis and kmeans clustering. Importantly: Rows must be named after features.
#' @param centrality_methods A character vector with strings specifying the method for selecting the most central feature of a cluster:
#' \itemize{
#'   \item \strong{two-in-a-row} - using PCA, selects the feature that shows up two times in a row as we calculate sum of squares adding more and more PC's is selected
#'   \item \strong{max-depth} - using PCA, selects the feature with the maximum sum of squares calculated across the number of pc's requested as the "max_depth"
#'   \item \strong{first-most-frequent} - using PCA, determines the max sum of squares for 2 pcs, 3 pcs, 4 pcs ... up to N pc's and then picks the feature that showed up the most times across all those calculations
#'   \item \strong{mhorn} - feature most similar to others (ie, largest sum to all other elements) wins
#'   \item \strong{spearman} - feature most similar to others (ie, largest sum to all other elements) wins
#'   \item \strong{pearson} - feature most similar to others (ie, largest sum to all other elements) wins
#'   \item \strong{by-rank} - defaults to the most significant according to \code{rank_df}
#' }
#' @param cluster_id_width An integer indicating how many characters to use for cluster group and cluster number id's. Defaults to one more than the number of characters in \code{max_clusters}.
#' @param cluster_plot_sizes Integer vector indicating which cluster groups to save as plots with clusters circled and central elements labeled. Only used if \code{centrality_methods} is one of the pca options.
#' @param file_prefix The text to be prepended to the file names for tables and plots
#' @param max_clusters Integer indicating the maximum number of clusters to split data into
#' @param min_clusters Integer indicating the minimum number of clusters to split data into
#' @param max_depth Integer indicating the maximum depth across principle components to use for determining most central element
#' @param my_core_number Integer value specifying to number of parallel processes to use when calculating mhorn indices. Defaults to 1.
#' @param my_seed The seed key to use so clustering can be reproduced
#' @param output_central_elements Boolean whether or not to save the table of central elements by cluster group
#' @param output_cumulative_variance Boolean whether to save a plot of the cumulative variance explained by the pca axes. Only used if \code{centrality_methods} is one of the pca options.
#' @param output_dir The base directory to which files and plots will be saved
#' @param output_gmt Boolean whether or not to save the gmt data to file
#' @param output_heatmap Boolean whether to save correlation heatmap to file. Ignored if \code{centrality_methods} is one of the PCA options.
#' @param output_pc1_vs_pv2 Boolean whether to save a plot of the principle component 1 and 2 axes. Only used if \code{centrality_methods} is one of the pca options.
#' @param output_ranked_central_elements Boolean whether to save the table of unique central elements sorted by rank within cluster group
#' @param rank_clm One-length character vector with the name of the column holding the initial rankings, if any, in either rank_df if one was sent, or in \code{feature_df} otherwise
#' @param rank_df Data.frame with \code{feature_df} features by row in column one and \code{rank_clm} with numeric default ranking for tie-breaking.  If \code{<NA>} \code{rank_clm} will be looked for in \code{feature_df}.
#' @param dist_method String indicating the method to pass to \code{stats::dist} method for clustering
#' @param hclust_method String indicating the method to pass to \code{stats::hclust} method for clustering
#'  
#' @return Returns 3 variable list with \code{cluster_members}, seed, and results. Results is a named list of each \code{centrality_methods} with \code{central_elements} and either pca or correlations ( depending on the \code{centrality_methods} )
#' 
#' @export
#' 
find_central_elements_by_cluster <- function( 
	feature_df, 
	centrality_methods ="max-depth", 
	cluster_id_width = NA,
	cluster_plot_sizes = NA, 
	file_prefix = "central_elements",
	max_clusters = NA,
	max_depth = NA, 
	min_clusters = 1L, 
	my_core_number = 1, 
	my_seed = NA, 
	output_central_elements = T,
	output_cumulative_variance = F,
	output_dir = ".",
	output_gmt = T,
	output_heatmap = F,
	output_pc1_vs_pc2 = F,
	output_ranked_central_elements = T,
	rank_clm ="Rank",
	rank_df = NULL,
	dist_method = "euclidean",
	hclust_method = "complete"
){
	
	rownames_have_numbers = suppressWarnings(any(!is.na(as.numeric(rownames(feature_df)))))
	if (rownames_have_numbers){
		stop("feature_df has numbers instead of feature names for rownames. Is feature_df a data.frame with features by row?")
	}
	
	if (is.na(cluster_id_width) ){
	  cluster_id_width <- ( as.character(max_clusters) %>% nchar )
	} 

	# ensure the rownames are set properly if we got a rank_df
	if ( !is.null(rank_df) ){
	  df_with_ranks <- rank_df
	  rownames(df_with_ranks) <- rank_df[[1]]
	} else {
		df_with_ranks <- feature_df
	}
	
	# get rank data from feature_df or rank_df using rank_clm
	rank_data <- get_ranks_from_df( df_with_ranks, rank_clm )
	# we must have rank values for all rows of feature_df
	if ( any( rownames(feature_df) %ni% names(rank_data)) ){
	  missing_ranks <- rownames(feature_df)[ which(rownames(feature_df) %ni% names(rank_data)) ]
	  warning( "The following features are missing from rank_df and will be dropped: ", paste0( missing_ranks, collapse=", " ))
	  rank_data = rank_data[names(rank_data) %ni% missing_ranks]
	}
	if( !is.na(rank_clm) ) feature_df[ rank_clm ] <- NULL #remove ranks from feature_df if they exist
	feature_df %<>% .[ names(rank_data)[names(rank_data) %in% rownames(feature_df)], ] #sort by rank_data order, 
	
	# drop the ranks that aren't used so the rank numbers don't blow out our plots
	rank_data = sort(rank_data[names(rank_data) %in% rownames(feature_df)])
	rank_data[] = 1:length(rank_data) 
	
	max_max_depth = nrow(feature_df)-1
	if (is.na(max_depth)){
		max_depth <- max_max_depth
	} else if ( max_depth > max_max_depth ){
		warning(paste0("max_depth was too large (", max_depth,")  and was set to one less than the total number of features (", max_max_depth, ")."))
		max_depth <- max_max_depth
	}
	
	if (is.na(max_clusters)){
		max_clusters <- max_max_depth
	} else if ( max_clusters > max_max_depth ){
		warning(paste0("max_clusters was too large (", max_clusters,")  and was set to one less than the total number of features (", max_max_depth, ")."))
		max_clusters <- max_max_depth
	}
	
	
	# remove non numeric columns
	operatable_clms = binfotron::operatable_columns(my_dt = data.table::as.data.table(feature_df), acceptable_classes = 'numeric')
	nonoperatable_clms = names(feature_df)[! names(feature_df) %in% operatable_clms]
	if (length(nonoperatable_clms) > 0){
		cat(paste0("Removing ", length(nonoperatable_clms), " column(s) that were not numeric: \n"))
		cat(paste0(nonoperatable_clms, collapse = '\n'))
		cat('\n')
		for (this_clm in nonoperatable_clms) feature_df[[this_clm]] = NULL
	}
	
	# go ahead and scale data since both pca and kmeans require it, 
	#   and it doesn't impact the various correlation methods
	# centering is false because mhorn can't handle negative values???
	# analysis2_df = scale(feature_df, center = F) %>% as.data.frame()
	# feature_df %<>% {as.data.frame(scale(., center = F))}
	# my_clm = unlist(analysis2_df[[1]])
	# > sqrt(sum(my_clm^2)/(length(my_clm)-1))
	# [1] 1
	cat("Scaling analysis data by feature (across samples).\n\n")
	feature_df %<>% t() %>% scale(center=F) %>% t() %>% data.frame()

	cat(paste0("Generating cluster groups from ", min_clusters, " to ", max_clusters, " clusters using seed ", my_seed, ".\n\n"))

	
	# wanted to be able to look at the results of the clusters alongside the values
	dd <- dist(feature_df, method = dist_method )
	my_hclust <- hclust(dd, method = hclust_method )

	should_continue = TRUE
	cluster_members = list()
	for (clust_index in 1:max_clusters){
		if ( !should_continue ) break;
		should_continue = tryCatch({
			cut_results = cutree(my_hclust, k = clust_index)
			for ( branch_index in 1: clust_index){
				branch_name =	paste0("Cluster_", stringr::str_pad(clust_index, width=cluster_id_width, pad="0")  , "_", stringr::str_pad(branch_index, width=cluster_id_width, pad="0"))
				my_names = names(cut_results[cut_results == branch_index])
				cluster_members[[branch_name]] = my_names
			}
			TRUE
		},
		error = function(cond){
			message("\nLast valid cluster was at at clust_index = ", clust_index - 1)
			FALSE
		})
	}

	rtn_list <- list(cluster_members=cluster_members, seed=my_seed, results=list())
	
	# ---------------------------------------------------------------------------
	# iterate over list of centrality_methods
	# ---------------------------------------------------------------------------
	for ( c_method in centrality_methods ) {
		#c_method =  'by-rank'
		
		cat(paste0(c_method, "\n"))
	  elements_data <- NA
  	if ( c_method == "mhorn" ){
  	  cat("Calculating Morisita-Horn similarity indices.\n\n")
  		if (any(feature_df < 0)) {
  			warning("MHorn cannot handle negative numbers so all number will be shifted up by the min value.")
  			feature_df = feature_df - min(feature_df)
  		}
  	  elements_data <- compare_via_mhorn(feature_df)
  	} else if (c_method %in% c("spearman", "pearson")){ 
  		#we have to transpose the input matrix because cor works on columns so we need elements as columns ...
  		cat(paste0("Calculating ", c_method, " correlation indices.\n\n"))
  		
  		elements_data <- cor(t(feature_df), use="pairwise.complete.obs", method=c_method)
  	} else if (c_method %in% c("by-rank")){ # by-rank will use the the correlation elements just for plotting
  		#we have to transpose the input matrix because cor works on columns so we need elements as columns ...
  		cat(paste0("Calculating correlation indices for plotting by-rank cluster results.\n\n"))
  		
  		elements_data <- cor(t(feature_df), use="pairwise.complete.obs", method="pearson")
  	} else {
  		cat("Running principle component analysis on centered data.\n\n")
  	  #do pca, data is already scaled but isn't centered
  	  pca <- prcomp(feature_df, scale=F, center = T)
  	  elements_data <- pca$x
  	}
	
  	# Find lowest ranked/most central element in each cluster
  	cat(paste0("Finding central elements in each cluster using ", c_method, " method\n\n"))
  
  	
  	central_elements_res <- isolate_central_cluster_elements( 
  		elements_data = elements_data, 
  		cluster_members = cluster_members, 
  		element_ranks = rank_data, 
  		max_depth = max_depth, 
  		centrality_method = c_method
  	)
  	central_elements = sapply(central_elements_res,function(x){x$central_element})
  	index_values = lapply(central_elements_res,function(x){x$all_values})
  	rm(central_elements_res)
  	
  	max_cluster_value = as.numeric(strsplit(names(index_values[length(index_values)]),split="_")[[1]][2])
  	
  	# put a table together of all the average index/correlation values for each feature
  	index_ave_mx = matrix(nrow = nrow(elements_data), ncol = max_cluster_value)
  	rownames(index_ave_mx) = rownames(elements_data)
  	colnames(index_ave_mx) = 	paste0("Cluster_", stringr::str_pad(1:max_cluster_value, width=cluster_id_width, pad="0"))
  	for (clm_name in colnames(index_ave_mx)){
  		my_branches = grep(paste0(clm_name,"_"),names(index_values), value = T)
  		for (my_branch in my_branches){
  			my_values = index_values[[my_branch]]
  			index_ave_mx[ names(my_values), clm_name ] = my_values
  		}
  	}
  	
  	#ensure output directory exists and create paths
  	full_dir <- file.path(output_dir, switch(c_method, `by-rank`="by_rank", `first-most-frequent`="pca_fmf", `max-depth`="pca_max", `two-in-a-row`="pca_tir", c_method ) )
  	dir.create(full_dir, showWarnings=F)
  	gmt_file_output_path = file.path(full_dir, paste0(file_prefix, "_gmt.txt"))
  	central_elements_output_path = file.path(full_dir, paste0(file_prefix, "_features.tsv"))
  	ranked_central_elements_output_path = file.path(full_dir, paste0(file_prefix, "_ranked.tsv"))
  	heatmap_output_path = file.path(full_dir, paste0(file_prefix, "_heatmap.pdf"))
  	heatmap_selection_output_path = file.path(full_dir, paste0(file_prefix, "_heatmap_selection.pdf"))
  	cumulative_variance_output_path = file.path(full_dir, paste0(file_prefix, "_variance.jpg"))
  	pc1_vs_pc2_output_path = file.path(full_dir, paste0(file_prefix, "_pc1_vs_pc2.jpg"))
  	
  	#write files if requested
    if ( output_gmt ){
      write_list_as_gmt( cluster_members, central_elements, gmt_file_output_path )
    }
    if ( output_central_elements ){
      write_central_elements_table( 
        central_elements, 
        rownames(feature_df), 
        output_path = central_elements_output_path, 
        cluster_id_width = cluster_id_width 
      )
    }
    if ( output_ranked_central_elements ){
      if ( all(is.na(rank_data)) ){
        warning( "Could not output the ranked central elements table because no rankings were provided as rank_clm or rank_df.")
      } else {
        write_ranked_central_elements_table( 
          central_elements, 
          rank_data, 
          output_path = ranked_central_elements_output_path 
        )
      }
    }
  	#do various plots
  	if (c_method %in% c("mhorn", "spearman", "pearson", 'by-rank') ){
  	  if ( output_heatmap ){
        method_formatted <- stringr::str_to_title(c_method)
        if (c_method %in% c("spearman", "pearson") ){
        	my_plot_title = paste(method_formatted, "Correlations")
        	legend_title = "Rho (abs)"
        } else {
        	my_plot_title = method_formatted
        	legend_title = "Similarity"
        }
        
        row_dend = as.dendrogram(my_hclust)
        hm_data = abs(elements_data)
        # make same feat comparisons NA's
        for (this_feature in rownames(hm_data)){
        	hm_data[this_feature, this_feature] = NA
        }
        
        grid_units = 'mm'
        grid_size = 100
        my_font_size = grid_size/nrow(hm_data) * 2.5
        
        col_fun <- circlize::colorRamp2(length(rank_data):1, plasma(length(rank_data)))
        row_ha = rowAnnotation(Rank = rank_data, col = list(Rank = col_fun), show_legend = FALSE)
        
        # drawing rect on complexheatmap: https://github.com/jokergoo/ComplexHeatmap/issues/522
        my_hm = Heatmap(
        	hm_data, 
        	name = 'hm_data',
        	col=viridis(100),
        	cluster_rows = row_dend, 
        	cluster_columns = row_dend,
        	show_column_dend = FALSE,
        	show_column_names = FALSE,
        	width = unit(grid_size, grid_units), 
        	height = unit(grid_size, grid_units),
        	column_names_gp = grid::gpar(fontsize = unit(my_font_size, grid_units)),
        	row_names_gp = grid::gpar(fontsize = unit(my_font_size, grid_units)),
        	row_dend_side = "right",
        	row_names_side = "left",
        	row_dend_width = unit(0.5, "cm"),
        	right_annotation = row_ha,
        	heatmap_legend_param = list(
        		title = legend_title, 
        		direction = "horizontal",
        	  legend_width = unit(grid_size*0.3, grid_units),
        	  legend_height = unit(grid_size*0.07, grid_units),
        	  title_position = "lefttop"
        	)
        )
        	
        pdf(file=heatmap_output_path, width = 11, height = 8)

	        my_hm = draw(my_hm, heatmap_legend_side = "bottom")
	        ro = row_order(my_hm)
	        co = column_order(my_hm)
	        row_names = rownames(hm_data)[ro]
	        
	        if (!is.null(cluster_members)){
	        	# figure out max number of clusters ( color will be based on that )
	        	max_clusters = max(sapply(names(cluster_members), function (x){as.numeric(strsplit(x, split = "_")[[1]][2])}) )
	        	my_colors = rev(viridisLite::plasma(max_clusters - 1)) # don't need to outline the whole plot or the diagonals
	        	max_nchar = nchar(max_clusters)
	        	tile_size = 1/length(ro)
	        	
	        	for (cluster_number in max_clusters:1){
	        		for (branch_number in 1:max_clusters){
	        			raw_branch_name = paste0("Cluster_", stringr::str_pad(cluster_number, width=cluster_id_width, pad="0")  , "_", stringr::str_pad(branch_number, width=cluster_id_width, pad="0"))
	        			branch_features = cluster_members[[raw_branch_name]]
	        			branch_features = row_names[row_names %in% branch_features]
	        			if ( length(branch_features) > 1) {
	        				# label the correlation raw values
		        			box_height = tile_size * length(branch_features) # same as width for raw values
		        			start_point = (which(row_names == branch_features[1])-1)/length(row_names)
		        			if (cluster_number > 1){ # no point in putting a box around everything
			        			decorate_heatmap_body("hm_data", row_slice = 1, column_slice = 1, {
			        				grid.rect(
			        					x = unit(start_point, "npc"),
			        					y = unit(1-start_point, "npc"), # top left
			        					width = box_height,
			        					height = box_height,
			        					gp = gpar(col = my_colors[cluster_number], fill = NA, lwd = 1, lty = 1), just = c("left", "top"))
			        			})
		        			}
		        			mean_branch_col_name = paste0("Branches ", stringr::str_pad(cluster_number, width=cluster_id_width, pad="0"))
	        			} # end if
	        			
	        		} # end for
	        	} # end for
	        	
	        } # if cluster_members
	        
        dev.off()
        
        pdf(file=heatmap_selection_output_path, width = 11, height = 8)

          
	        colnames(index_ave_mx) = gsub("Cluster_", "Branches ", colnames(index_ave_mx))
	        my_values_hm = Heatmap(
	        	index_ave_mx[ro,],
	        	name='index_ave',
	        	col=viridis(100),
	        	show_heatmap_legend = FALSE,
	        	width = unit(grid_size, grid_units), 
	        	height = unit(grid_size, grid_units),
	        	cluster_rows = FALSE,
	        	cluster_columns = FALSE,
	        	column_names_gp = grid::gpar(fontsize = unit(my_font_size, grid_units)),
	        	row_names_gp = grid::gpar(fontsize = unit(my_font_size, grid_units))
	        )
	        my_values_hm = draw(my_values_hm, heatmap_legend_side = "bottom")
	        
	        # box around the branch id's 
	        # mark the best with a red dot
	        
	        
	        if (!is.null(cluster_members)){
	        	# figure out max number of clusters ( color will be based on that )
	        	max_clusters = max(sapply(names(cluster_members), function (x){as.numeric(strsplit(x, split = "_")[[1]][2])}) )
	        	my_colors = rev(viridisLite::plasma(max_clusters - 1)) # don't need to outline the whole plot or the diagonals
	        	max_nchar = nchar(max_clusters)
	        	tile_size = 1/length(ro)
	        	
	        	for (cluster_number in max_clusters:1){
	        		for (branch_number in 1:max_clusters){
	        			raw_branch_name = paste0("Cluster_", stringr::str_pad(cluster_number, width=cluster_id_width, pad="0")  , "_", stringr::str_pad(branch_number, width=cluster_id_width, pad="0"))
	        			branch_features = cluster_members[[raw_branch_name]]
	        			branch_features = row_names[row_names %in% branch_features]
	        			if ( length(branch_features) > 1) {
	        				# label the correlation raw values
	        				box_height = tile_size * length(branch_features) # same as width for raw values
	        				start_point = (which(row_names == branch_features[1])-1)/length(row_names)

	        				mean_branch_col_name = paste0("Branches ", stringr::str_pad(cluster_number, width=cluster_id_width, pad="0"))
	        				
	        				if( mean_branch_col_name %in% colnames(index_ave_mx) ){
	        					index_ave_width = 1/ncol(index_ave_mx)
	        					
	        					# need the index of the best feature
	        					central_element = central_elements[raw_branch_name]
	        					text_start_point = (which(row_names == central_element))/length(row_names)
	        					
	        					
	        					decorate_heatmap_body("index_ave", row_slice = 1, column_slice = 1, {
	        						grid.rect(
	        							x = unit((cluster_number-1)/ncol(index_ave_mx), "npc"),
	        							y = unit(1-start_point, "npc"),
	        							width = unit(index_ave_width, "npc"),
	        							height = box_height,
	        							gp = gpar(col = "white", fill = NA, lwd = 1, lty = 1), just = c("left", "top")
	        						)
	        						grid.text(
	        							"*", 
	        							gp = gpar(col = "red"),
	        							x = unit((cluster_number - 0.5)/ncol(index_ave_mx), "npc"), 
	        							y = unit(1-text_start_point, "npc"), 
	        							rot = 0
	        						)
	        					})
	        				}
	        			} # end if
	        			
	        		} # end for
	        	} # end for
	        	
	        } # if cluster_members
        
        dev.off()
        
      } # if output_heatmp
  	} else if ( c_method != "by-rank" ) { #pca plots
  		
      #plot cumulative variance by principle component if requested
      if ( output_cumulative_variance){
  
        var_pt = pca_variance_plot(
        	pca_results = pca, 
        	output_path = cumulative_variance_output_path
        	#significance_threshold = cumulative_variance_threshold
        )
      }
      
      #plot pc1 vs pc2 as scatterplot if requested
      if ( output_pc1_vs_pc2 ) {
      	pc1v2_pt = pc1_vs_pc2_plot( elements_data %>% as.data.frame(), pc1_vs_pc2_output_path )
      }
      #plot clusters as individual plots for each # clusters requested
      if ( !all(is.na(cluster_plot_sizes) ) ){
        cat("Generating cluster plots for cluster groups ", cluster_plot_sizes, ".")
        #reformat the long list of cluster_members as data.frame with all elements ( i.e. features ) as rows and a pair of columns per cluster count indicating cluster membership and central elements 
        clusters_by_ct_df <- cluster_list_to_df( cluster_members, central_elements )
        elements_df <- elements_data %>% as.data.frame()
        #loop over various numbers of clusters requested to be plotted and plot them
        cluster_grps_to_plot <- paste("Clusters", stringr::str_pad(cluster_plot_sizes, width=cluster_id_width, pad="0"), sep="_") 
        for ( cluster_grp in cluster_grps_to_plot ){
          cluster_group_plot( 
            elements_df, 
            clusters_by_ct_df[[cluster_grp]], 
            central_elements_by_cluster = clusters_by_ct_df[[paste(cluster_grp, "cfs", sep="_")]], 
            cluster_group_id = cluster_grp, 
            output_path = file.path(full_dir, paste0(file_prefix, "_", tolower(cluster_grp), "_seed", my_seed, ".pdf")) 
          )
        }
      }
	  }
  	if ( c_method %in% c("mhorn", "spearman", "pearson") ){
  	  rtn_list$results[[ c_method ]] <- list( correlations=elements_data %>% as.data.frame(), central_elements=central_elements )
  	} else if (c_method != "by-rank" ){ 
  	  rtn_list$results[[ gsub( "-", "_", c_method ) ]] <- list( pca=pca, central_elements=central_elements )
  	} else{
  	  rtn_list$results[[ gsub( "-", "_", c_method ) ]] <- list( rankings=rank_data, central_elements=central_elements )
  	}
	}
	
	return(rtn_list)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate_kmeans_cluster_list
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Obtain a list of clusters with their elements from data matrix sent
#' 
#' @description Runs kmeans method requesting each number of clusters from min_clusters:max_clusters
#' 
#' 
#' @param elements_df Data.frame or matrix of data to be clustered as samples x elements ( i.e. features )
#' @param cluster_id_width The number of characters to include in cluster group and individual cluster id's ( will be used to left-pad cluster numbers with leading 0's )
#' @param max_clusters Integer number representing largest number of clusters to split data into. Must be less than number of rows in matrix.
#' @param min_clusters Integer number representing smallest number of clusters to split data into
#' @param my_seed The seed key to use before each call to kmeans method so each run can be reproduced
#'  
#' @return Returns lists of clusters named as cluster_#clusters-in-run_cluster#-from-run = c(cluster_variable_names). Element names in clusters are sorted alphabetically.
#' 
#' @export
#'
generate_kmeans_cluster_list = function(
	elements_df, 
	cluster_id_width = NA,
	max_clusters = nrow(elements_df)-1, 
	min_clusters = 2, 
	my_seed = NA
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# compare_via_mhorn
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Compare rows pairwise using Morisita Horn similarity index
#' 
#' @description Compare each pair of rows in the feature_df using Morisita Horn method and return a correlation matrix ( as data.frame )
#' 
#' 
#' @param feature_df Data.frame with rows to be compared pair-wise
#' @param my_core_number Integer number representing the number of parallel processes to use for mhorn calculations
#'  
#' @return Returns a correlation matrix in data.frame format
#' 
#' @export
#'
compare_via_mhorn = function( 
  feature_df, 
  my_core_number=4 
){
  #todo: potentially save mhorn analysis data to file and allow reloading it instead of reprocessing ... would need to add unique save path with seed ...
  other_elements <- rownames(feature_df)
  num_elements <- nrow(feature_df)
  my_core_number %<>% min(nrow(feature_df))
  all_results <- parallel::mclapply(rownames(feature_df)[-nrow(feature_df)], function(element_a){
    other_elements <<- other_elements[-1]
    correlations <- mapply(function(element_b){
      print(paste(element_a, ":", element_b))
      sub_feature_df <- feature_df[c(element_a, element_b),]
      my_empty_cols <- sapply(sub_feature_df, function (k) all(is.na(k)))
      sub_feature_df = sub_feature_df[, !my_empty_cols]
      my_result = NA
      if(ncol(sub_feature_df)>0){
        sub_feature_df[is.na(sub_feature_df)] <- 0 # make this mclapply
        if(all(rowSums(sub_feature_df, na.rm=T)>0)) {
          my_result <- vegan::vegdist(sub_feature_df, method="horn", binary=FALSE, diag=FALSE, upper=TRUE, na.rm = FALSE) %>% as.matrix()
          my_result = 1-my_result[1,2] #convert to similarity index ...
        } else if(any(rowSums(sub_feature_df, na.rm=T)==0)){ # if one of these is zero there is no similarity
          my_result = 0
        }
      }
      print(my_result)
      my_result
    }, other_elements)
    return(correlations)
  },mc.cores=my_core_number)
  
  #now populate the symmetrical dataframe to return
  mhorn_df <- data.frame(matrix(0,nrow=nrow(feature_df), ncol=nrow(feature_df)))
  names(mhorn_df) <- rownames(feature_df)
  rownames(mhorn_df) <- rownames(feature_df)
  for (f_i in 1:length(all_results)){
    corrs <- all_results[[f_i]]
    rep_indices <- (num_elements-length(corrs)+1):num_elements
    mhorn_df[rep_indices, f_i] <- corrs
    mhorn_df[f_i,rep_indices] <- corrs
  }
  mhorn_df
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# isolate_central_cluster_elements
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Determine most representative element in each cluster based on data provided
#' 
#' @description Takes pca scores or a correlation matrix ( as data.frame ) and uses it to determine the most representative element from each of a list of clusters. Clusters with 2 elements use ranks sent or random selection to determine centrality while clusters larger than 3 use the centrality_method
#'
#' 
#' @param elements_data A data.frame or matrix containing either elements x principle components ( as scores/x from PCA Analysis ) OR similarity scores of elements x elements
#' @param cluster_members A named list of clusters with their elements
#' @param element_ranks A named integer vector indicating the initial element rankings to be used for selection of best elements in clusters of length 2
#' @param max_depth An integer indicating the maximum number of Principle Components to use in determining best elements. Only used for PCA type centrality_methods
#' @param centrality_method A character vector with strings specifying the method for selecting the most central feature of a cluster:
#' \itemize{
#'   \item \strong{two-in-a-row} - using PCA, selects the feature that shows up two times in a row as we calculate sum of squares adding more and more PC's is selected
#'   \item \strong{max-depth} - using PCA, selects the feature with the maximum sum of squares calculated across the number of pc's requested as the "max_depth"
#'   \item \strong{first-most-frequent} - using PCA, determines the max sum of squares for 2 pcs, 3 pcs, 4 pcs ... up to N pc's and then picks the feature that showed up the most times across all those calculations
#'   \item \strong{mhorn} - feature most similar to others (ie, largest sum to all other elements) wins
#'   \item \strong{spearman} - feature most similar to others (ie, largest sum to all other elements) wins
#'   \item \strong{pearson} - feature most similar to others (ie, largest sum to all other elements) wins
#'   \item \strong{by-rank} - defaults to the most significant according to \code{rank_df}
#' }
#' @return Returns original list of cluster_members with the most representative element for each named cluster
#' 
#' @export
#'
isolate_central_cluster_elements = function(
	elements_data, 
	cluster_members, 
	element_ranks=NA,
	max_depth=NA, 
	centrality_method="max-depth" 
){
	
	if (is.na(max_depth) ){
		max_depth = ncol(elements_data)-1
	}

  #convert incoming data.frame to matrix for speed optimization
  if ( is.data.frame(elements_data) ){
    clm_names <- colnames(elements_data)
    elements_data %<>% as.matrix()
    colnames(elements_data) <- clm_names
    rm(clm_names)
  }
	max_depth %<>% min( (ncol(elements_data) - 1), na.rm = T )
  has_ranks = T
	if ( all(is.na(element_ranks)) ){
	  cat("element_ranks was not provided. Note that this means clusters with 2 elements will select the highest ranked element.\n\n")
	  has_ranks = F
	}
	#iterate over clusters in cluster list and for each cluster return the most central
	central_element_list<-lapply(cluster_members, function(cluster){
		# print(cluster)
		all_values = NA
		#x should be a character vector of length = #elements in this cluster
		#determine the most representative element
		if ( length(cluster) == 1 ){
			central_element <- cluster[1] 
			#for clusters of 1 element, return the element
		} else if ( length(cluster) == 2 ){
			central_element <- cluster
		} else if ( length(cluster) > 2 ){
			#for clusters over 2, return the most central element via centrality_method requested
		  if ( centrality_method %in% c("mhorn", "spearman", "pearson", 'by-rank') ){
		  	# for by-rank we only need these data to display correlations on the heatmaps
		    central_element_res <- correlation_central_element( 
		    	correlation_data = elements_data, 
		    	elements = cluster 
		    )
		    if ( centrality_method == "by-rank"){
		    	central_element <- element_ranks[ cluster ] %>% { names(.)[ which.min(.) ] }
		    }	else {
		    	central_element = central_element_res$most_similar
		    }
		    all_values = central_element_res$values
		  } else {
			  central_element <- pca_central_element( 
			  	pca_data = elements_data,
			  	elements = cluster,
			  	max_depth = max_depth,
			  	centrality_method = centrality_method 
			  )
		  }
		}
		
		# correlation_central_element will return ties so we need to break those
		if ( length(central_element) > 1 ){ # for clusters of 2 elements, return the lowest ranked element OR random selection if no element ranking was sent
			if ( has_ranks ){
				if ( element_ranks[[central_element[1]]] < element_ranks[[central_element[2]]] ) {
					central_element <- central_element[1]
				} else {
					central_element <- central_element[2]
				}
			} else {
				#return first element in cluster
				central_element <- central_element[1]
			}
		}
		
		
		return(list(central_element = central_element, all_values = all_values))
	})
	#print(paste("Total time:", Sys.time()-time_in))
	return(central_element_list)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pca_central_element
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Determine the most central of a list of elements based on pca scores
#' 
#' @description Calculates maximum total euclidean distance from each element to all other elements through as many principle components as necessary based on the method indicated.
#' 
#' 
#' @param pca_data Data.frame or matrix of pca loadings ( elements x principle components )
#' @param elements Character vector of the set of elements from which to determine centrality. Values must be a subset of pca_data rownames.
#' @param max_depth An integer indicating the maximum number of Principle Components to use in determining most central element
#' @param centrality_method A one-length character vector with value "two-in-a-row", "max-depth" or "first-most-frequent". "two-in-a-row" selects whichever feature shows up two times in a row as you compare more and more pc's. So if you get feature a when including only pc1-pc2 , feature b for pc1-pc3, and feature b for pc1-pc4, it would return feature b.  "first-most-frequent" returns whichever feature shows up the most times as you compare each number of PC's. "First" because if different features each show up the same number of times whichever showed up with the fewest pc's is returned. 
#'  
#' @return Returns the element name that is most central based on minimum euclidean distance through n principle components.
#' 
pca_central_element <- function( 
	pca_data, 
	elements, 
	max_depth = ncol(pca_data)-1, 
	centrality_method = "max-depth" 
){
  if( !centrality_method %in% c("max-depth", "first-most-frequent", "two-in-a-row") ){
    stop(paste("Centrality method '", centrality_method, "' not recognized. Central element can not be obtained."))
  }
  #number of principle component vectors to include in "most central" calculation
  if ( max_depth < 2 ) stop("The max_depth must be 2 or greater.")
  #max_depth can't be greater than the number of PC's
  max_depth %<>% min( (ncol(pca_data)-1), na.rm=T )

  # convert to matrix if this is a df - matrix is DRAMATICALLY faster for indexing
  if ( is.data.frame(pca_data) ){
    clm_names <- names(pca_data)
    pca_data %<>% as.matrix()
    colnames(pca_data) <- clm_names
  }
	#init vectors to hold the calculated euclidean distances and central elements
	euclidean_distances <- vector(mode = "numeric", length = length(elements))
	names(euclidean_distances) <- elements
	all_central_elements = vector(mode = "character", length = max_depth)
	if ( centrality_method == "max-depth" ) {
		depth <- max_depth
	} else {
		depth <- 2
	}
	while ( depth <= max_depth ){
		euclidean_distances[euclidean_distances != 0] <- 0
		pcs <- c(1:depth)
		for ( element in elements ){
			#calculate distance to each of the other elements
			tot_distance <- 0
			e1_pcs <- pca_data[ rownames(pca_data) == element, ] %>% unlist()
			for ( other_element in elements[-which(elements == element)] ){
				sum_of_squares <- 0
				e2_pcs <- pca_data[ rownames(pca_data) == other_element, ] %>% unlist()
				#along k principle components
				for ( pc_axis in pcs ){
					sum_of_squares <- sum_of_squares + (e2_pcs[ pc_axis ] - e1_pcs[ pc_axis ])^2
				}
				tot_distance <- tot_distance + sqrt(sum_of_squares)
			}
			euclidean_distances[element]<-tot_distance
		}
		all_central_elements[depth] <- names(which.min(euclidean_distances)) #names(euclidean_distances)[which.min(euclidean_distances)]
		if ( centrality_method == "two-in-a-row" && all_central_elements[depth-1] == all_central_elements[depth] ) break #exit if we have found the same "central_element" twice in a row as we increase depth
		depth <- depth + 1
	}
#	print(paste("Length ", length(elements), "cluster"))
	if ( centrality_method == "first-most-frequent" ){
	  return(sort(table(all_central_elements), decreasing=T) %>% {attr(., "dimnames")$all_central_elements[1]})
	} else {
		return(all_central_elements[min(depth,max_depth)])
	}
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# correlation_central_element
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Determine the most similar of a list of elements to all other elements in the list based on a correlation matrix
#' 
#' @description Calculates maximum similarity/correlation between any single element of group and all other elements
#' 
#' 
#' @param correlation_data Data.frame or matrix of similarity scores / correlations containing at a minimum all elements. This df should be symmetrical - i.e. have equivalent values for row A to column B and row B to column A
#' @param elements Character vector of the set of elements ( i.e. features ) from which to determine greatest similarity. Values must be a subset of correlation_data row and column names
#'  
#' @return Returns the name of the element that is most similar to all other elements in the set
#' 
correlation_central_element <- function(
  correlation_data, 
  elements
){
  if (is.data.frame(correlation_data) ){
    clm_names <- names(correlation_data)
    correlation_data %<>% as.matrix()
    colnames(correlation_data) <- clm_names
  }
  #sanity check that all rownames exist as colnames in correlation_data
  if (!all(colnames(correlation_data) %in% rownames(correlation_data)) | !all(elements %in% colnames(correlation_data))){
    warning("All elements must exist in both colnames and rownames of correlation_data - returning NULL.")
    return(NULL)
  }
  #sanity check whether correlation_data is symmetrical with regard to similarity values
  element_1 <- colnames(correlation_data)[1]
  element_2 <- colnames(correlation_data)[2]
  if (!identical(correlation_data[element_1, element_2], correlation_data[element_2, element_1]) | is.na(correlation_data[element_1, element_2])){
    warning("Values for correlation_data[a,b] must equal values for correlation_data[b,a], and should not be NA.")
  }

  #we have similarity indices ( in correlation_data ), so we just need to average these similarity indices from each element to all other elements
  #the element with the largest mean to all other elements is the "most similar" element of the cluster
  
  # SKIP NA's in mean
  # SKIP comparison with self in mean
  my_means = c()
  for (element in elements){
  	mean_abs_correlation = mean(abs(unlist(correlation_data[element, elements[elements != element]])), na.rm = T)
  	my_means = c(my_means, mean_abs_correlation)
  }
  names(my_means) = elements
  most_similar = my_means[my_means == max(my_means)] # we'll handle duplicates later
  
  return(list(most_similar = names(most_similar), values = my_means))

}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# cluster_list_to_df
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Reformat a long list of cluster members as a data.frame
#' 
#' @description Reformats long-form list of cluster members as a data.frame with all elements as rows and one or two columns per cluster indicating cluster membership and optionally central elements for each cluster 
#' 
#' 
#' @param cluster_members Long list of each cluster in each cluster count and the elements that belong to that cluster ( i.e. Cluster_002_001 <- Bindea_aDC, Bindea_BCells ... ). 
#' @param central_elements Optional list of central elements by cluster id ( i.e. Cluster_002_001 <- central element of that cluster ). If omitted, central element columns will not be added to return data.frame.
#'  
#' @return Returns the data.frame containing cluster membership for each group of clusters in the cluster_members list and optionally the central element of each cluster.
#' 
cluster_list_to_df <- function( 
  cluster_members, 
  central_elements=NA 
){
  all_elements <- unique(unlist(cluster_members)) # find full list of elements
  return_df <- data.frame( row.names=all_elements )
  for ( cluster_key in names(cluster_members)){
    #get number of clusters and cluster number
    cluster_group_id <- paste("Clusters", strsplit(cluster_key,split = "_")[[1]][2], sep="_")
    if( !cluster_group_id %in% names(return_df) ){
      return_df[cluster_group_id] <- 0
      if( !all(is.na(central_elements) ) ){
        clusters_cf_id <- paste(cluster_group_id, "cfs", sep="_")
        return_df[ clusters_cf_id  ] <- ""
      }
    } 
    cluster_num <- strsplit(cluster_key, split = "_")[[1]][3] %>% as.integer()
    for ( element in cluster_members[[ cluster_key ]] ){
      return_df[ element, cluster_group_id ] <- cluster_num
      if ( !all(is.na(central_elements)) && central_elements[ cluster_key ] == element ){
        return_df[ element, clusters_cf_id ] <- element
      }
    }
  }
  # make all non central element columns factors
  non_cfs_clms <- names(return_df)[!names(return_df) %like% "_cfs"]
  for (cluster_col in non_cfs_clms){
    return_df[[ cluster_col ]] %<>% as.factor()
  }

  return(return_df)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write_list_as_gmt
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Saves to file list of element combinations and most central element from each grouping
#' 
#' @description Creates a file with one tab separated line per each cluster of elements formatted as central_element__clusterid, comment, full list of elements in this gene signature with NA's for extra spots to length of largest cluster
#' 
#' 
#' @param cluster_members Named list of clusters ( names will be used as clusterid's )
#' @param central_elements Named list of central elements ( i.e. features ) from each cluster ( same length as cluster_members with matching names )
#' @param output_path One-length character vector with fully qualified path to output file
#' @param empty_cell Filler for empty cells. Typically '' or NA. Note there is a gene called 'NA'.  It's in the influenza genome but still...
#' @param display_output Boolean option to display which elements are in the final signatures.  For gene signatures this will typically be too much information to display
#'  
#' @return Returns gene signature data as saved to file
#' 
write_list_as_gmt = function( 
	cluster_members, 
	central_elements, 
	output_path = file.path(getwd(),"output_file.gmt.txt"),
	empty_cell = '',
	display_output = F
){
	#create output as concatenated central_element__clusterid, comment, full list of elements in this gene signature with empty_cell's for extra spots
	max_cluster_len <- max(mapply(length, cluster_members))
	gene_sig_output <- c()
	for (cluster in seq_along(cluster_members)){
		cluster_name <- names(cluster_members)[cluster]
		extras_to_add <- c(rep(empty_cell,max_cluster_len - length(cluster_members[[cluster_name]])))
		gene_sig_output %<>% c(paste(paste(cluster_name, central_elements[[cluster_name]], sep="__"), 
																	 paste("length", length(cluster_members[[cluster_name]]), "cluster"), 
																	 paste0(c(cluster_members[[cluster_name]], extras_to_add), collapse="\t"),sep="\t"))
	}
	
	gene_sig_output %<>% paste0(collapse="\n") 
	writeLines(gene_sig_output, output_path)
	
	cat(paste0("\nGene signature data saved to:", output_path, '\n\n'))
	if (display_output) gene_sig_output
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write_central_elements_table
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Saves to file a dataframe of elements ( i.e. features ) x cluster groups
#' 
#' @description Creates a file with one tab separated row per each element and one column for each number of clusters with 1's for the central elements of each cluster in each cluster group and "" everywhere else
#' 
#' 
#' @param all_elements Character vector of element names containing, at a minimum, all central elements defined in first argument
#' @param central_elements Named list of central element from each cluster ( as returned by isolate_central_cluster_elements )
#' @param cluster_id_width The number of characters in cluster group id's and individual cluster id's
#' @param output_path One-length character vector with fully qualified path to output file
#' @param return_output Boolean option to display final output
#'  
#' @return Returns central elements data.frame as saved to file
#' 
write_central_elements_table = function(
	central_elements, 
	all_elements, 
	output_path=file.path(getwd(), "central_elements_by_num_clusters.tsv"),
	return_output = F,
	cluster_id_width = NA
){
  #get max and min number of clusters
  cluster_cts <- sort(names(central_elements))[c(1,length(central_elements))] %>% strsplit("_")
  min_clusters <- as.integer(cluster_cts[[1]][2])
  max_clusters <- as.integer(cluster_cts[[2]][2])

  cluster_id_width %<>% max(., nchar(cluster_cts[[1]][[2]]), na.rm=T)
  
  #create empty dataframe with correct number of rows and columns, then set colnames and rownames to cluster_id's and elements
  output_df <- data.frame(matrix("",ncol=max_clusters-min_clusters+1, nrow=length(all_elements)))
  names(output_df) <- paste("Clusters", stringr::str_pad(min_clusters:max_clusters,width = cluster_id_width, pad="0"), sep="_")
  rownames(output_df) <- all_elements
  
  for ( cluster in stringr::str_pad(min_clusters:max_clusters, width=cluster_id_width, pad="0") ){
    cluster_indices <- paste("Cluster", cluster, sep="_") %>% grep(names(central_elements))
    if( length(cluster_indices) == 0 ) next
    #pull the central elements for all clusters within this cluster group
    cluster_elements <- central_elements[ cluster_indices ] %>% sapply(function(x) return(x[[1]]) )
    #put 1's at the intersection of the cluster_elements and the cluster group which is, for instance, Clusters_002 rather than Cluster_002 (note the plural)
    output_df[ cluster_elements, names(output_df) == paste("Clusters", cluster, sep="_") ] <- 1
  }
  
  output_df = data.frame(Feature = rownames(output_df), output_df)
  rownames(output_df) = NULL
  
  if ( !is.na(output_path) ){
  	fwrite(output_df, file=output_path, sep = "\t")
    cat(paste0("Central elements table saved to:", output_path,'\n\n'))
  }
  
  if (return_output) return(output_df)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write_ranked_central_elements_table
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Saves to file a dataframe with two columns, Element and Rank
#' 
#' @description Creates a file with one tab separated row per each unique central cluster element along with its ranking
#' 
#' 
#' @param central_elements Named list of central element from each cluster ( as returned by isolate_central_cluster_elements ), assumes naming convention of {cluster_group}_{cluster_number_within_group}
#' @param element_ranks Named integer vector of element rankings
#' @param output_path One-length character vector with fully qualified path to output file
#' @param display_output Boolean option to display final output
#'  
#' @return Returns data.frame of unique ranked central cluster elements as saved to file
#' 
write_ranked_central_elements_table = function(
    central_elements, 
    element_ranks, 
    output_path=file.path(getwd(), "ranked_unique_central_elements_ranked.tsv"),
    display_output=F
){
  #find the first cluster that contains each unique central element
  unique_central_elements <- central_elements %>% unlist() %>% {.[!duplicated(.)]}
  #get the cluster_ids for the unique central elements found above after stripping the _{cluster_number_within_group}
  cluster_grp_ids <- names(unique_central_elements) %>% gsub("_[0-9]+$", "", .) %>% unique()
  ordered_central_elements <- c()
  #loop over the cluster id's that contain unique elements
  for ( cluster_grp in cluster_grp_ids ){
    #get the unique elements that are from clusters in this cluster group
    elements_by_cluster <- unique_central_elements[ names(unique_central_elements) %like% cluster_grp ]
    #because there could be multiple clusters in this cluster group with unique central elements, we need to order them by rank
    ordered_central_elements %<>% c( element_ranks[ elements_by_cluster ] %>% sort() %>% names() )
  }
  # create the output data frame with the ordered central elements and their ranks
  output_df <- data.frame(element=ordered_central_elements, rank=element_ranks[ordered_central_elements])
  write.table( x = output_df, file= output_path, sep="\t", row.names = F, col.names=T )
  cat(paste0("Ranked central elements table saved to:", output_path,'\n\n'))
  
  if (display_output) output_df
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pca_variance_plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Saves to file a cumulative variance plot for principle components of a PCA analysis
#' 
#' @description Saves to file a cumulative variance plot for principle components of a PCA analysis with optional coloring by cumulative percent variance explained
#' 
#' 
#' @param pca_results A prcomp object
#' @param output_path One-length character vector with fully qualified path to output file
#' @param my_width One-length numeric vector for the width of the output plot
#' @param my_height One-length numeric vector for the height of the output plot
#'  
#' @return Returns the ggplot object
#' 
#' @export
pca_variance_plot <- function( 
	pca_results, 
	output_path=NA, 
	my_width=8, 
	my_height=8 
){
  eigenvalues <- pca_results$sdev^2
  plot_df <- data.frame( pc=1:length(eigenvalues), percent=eigenvalues/sum(eigenvalues)*100)
  plot_df$cumulative_variance <- cumsum(plot_df$percent)
  #store number of principle components here for repeated use following
  num_pcs <- length(eigenvalues) 
  #depending how many principle components there are, set how many breaks ( and thus, labels ) to show
  if ( num_pcs <= 40 ) {
  	x_breaks <- 1:num_pcs
  } else {
  	x_breaks <- seq(0, num_pcs, round(num_pcs/40))
  }
  #plot it ...
  v_plot <- ggplot(plot_df, aes(factor(pc), percent))+#, color=thresh, fill=thresh) + 
    geom_bar(stat="identity") + 
    geom_path(aes(y=cumulative_variance), group=1) + 
    geom_point(aes(y=cumulative_variance)) + 
    scale_x_discrete(breaks=x_breaks) + 
    scale_y_continuous(n.breaks=10, breaks=waiver()) + 
    labs(x="Principle Component", 
         y="Percent of Variance Explained", 
         title="Percent Variance by Principle Component", 
         color="Cumulative\nVariance", 
         fill="Cumulative\nVariance") + 
    theme(plot.title=element_text(hjust=0.5))
  
  if( !is.na(output_path) ){
    ggsave(output_path, v_plot, width = my_width, height = my_height)
  	cat(paste0( "Variance plot saved to: ", output_path, ".\n\n"))
  }
  
  return(v_plot)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# cluster_group_plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Saves to file a plot of the clusters within a cluster group
#' 
#' @description Saves to file a plot showing the clusters within a cluster group, colored by group membership and with the central elements labeled for each cluster
#' 
#' 
#' @param scores_df Dataframe with elements ( i.e. features ) x principle component scores. Expects at least a PC1 and PC2 column
#' @param cluster_membership Integer factor with cluster membership for each element in scores_df
#' @param central_elements_by_cluster Character vector with element names at locations of most central elements for each cluster and "" everywhere else
#' @param cluster_group_id One-length character vector name for this cluster group to be used in plot title
#' @param output_path Character vector with fully qualified path to which plot will be saved
#' @param my_width One-length numeric vector with width of output plot ( in inches )
#' @param my_height One-length numeric vector with height of output plot ( in inches )
#'  
#' @return Returns the ggplot object
#' 
#' @export
cluster_group_plot <- function( 
  scores_df, 
  cluster_membership, 
  central_elements_by_cluster, 
  cluster_group_id, 
  output_path=NA, 
  my_width=8, 
  my_height=8 
){
  cluster_plot <- ggplot(scores_df, aes(PC1, PC2, color=cluster_membership))+
    geom_point() +
    labs(title=paste("Central elements of ", cluster_group_id), color="Cluster group") +
    ggrepel::geom_label_repel(aes(label=central_elements_by_cluster), size=1.5, min.segment.length = 0) +
    ggforce::geom_mark_ellipse() + #stat_ellipse()
    theme(plot.title=element_text(hjust=0.5))
  
  if(!is.na(output_path)){
    ggsave(output_path, cluster_plot, width = my_width, height = my_height)
  	cat(paste0( "Cluster group plot for ", cluster_group_id, " saved to: ", output_path, ".\n\n"))
  }
  
  return(cluster_plot)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pc1_vs_pc2_plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Saves to file a plot of pc1 vs pc2 from a pca analysis
#' 
#' @description Saves to file a scatterplot showing PC1 vs PC2
#' 
#' 
#' @param scores_df Dataframe with elements x principle components. Expects at least a PC1 and PC2 column
#' @param output_path Character vector with fully qualified path to which plot will be saved
#' @param my_width One-length numeric vector with width of output plot ( in inches )
#' @param my_height One-length numeric vector with height of output plot ( in inches )
#'  
#' @return Returns the ggplot object
#' 
pc1_vs_pc2_plot <- function( 
  scores_df, 
  output_path=NA, 
  my_width=8, 
  my_height=8 
){
  #scatterplot of scores for PC1 vs PC2
  pc_plot <- ggplot(scores_df,aes(PC1, PC2))+geom_point()+
    ggrepel::geom_label_repel(aes(label=rownames(scores_df)), size=1.25)+
    labs(title="Principle Components 1 vs. Principle Component 2") +
    theme(plot.title=element_text(hjust=0.5))
  
  if( !is.na(output_path) ){
    ggsave(output_path, width=my_width, height=my_height)
  	cat(paste0( "PC1 vs. PC2 plot saved to: ", output_path, ".\n\n"))
  }
  
  pc_plot
}

  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# group_clusters_by_size
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Group a list of clusters by size
#' 
#' @description Takes a list of clusters and groups them into named list by size
#' 
#' 
#' @param c_list Named list of clusters - names will be preserved in groupings
#' @param min_size Integer number representing smallest cluster size to keep
#' @param max_size Integer number representing largest cluster size to keep
#' @param include_empty_groups Boolean whether or not to include cluster sizes in return value for which there are no representatives
#'  
#' @return Returns lists of clusters ( with their elements ) grouped by cluster size named as clusters_of_x where x = cluster size.
#' 
#' @export
#'
group_clusters_by_size <- function(
  c_list, 
  min_size=1L, 
  max_size=NA, 
  include_empty_groups=F
){
  if ( !is.integer(min_size) ) stop("You must send an integer min_size.")
  #get vector of cluster sizes found in c_list
  sizes <- sort(unique(mapply(length, c_list)))
  #find the maximum cluster size if no max_size was sent
  if (is.na(max_size)) max_size = max(sizes)
  #subset sizes to be within min_size and max_size
  sizes %<>% .[ . >= min_size & . <= max_size ]
  #initialize list and create names
  clusters_by_size <- lapply(sizes, function(x) return( list() ) )
  names(clusters_by_size) <- c(paste("clusters_of", sizes, sep="_"))
  #subset c_list to only include clusters within the min_size:max_size criteria
  c_list %<>% .[ mapply(function(x) length(x) %in% sizes, .)]
  #loop over all clusters and place them in appropriate cluster_size list element
  for (c_index in 1:length(c_list)){
    c_size = length(c_list[[c_index]])
    clusters_by_size[[paste0("clusters_of_",c_size)]][[names(c_list)[c_index]]] <- c_list[[c_index]]
  }
  if ( !include_empty_groups ){
    clusters_by_size %<>% .[sapply(., function(c) length(c) > 0)]
  }
  return(clusters_by_size)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# collapse_identical_clusters
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Remove identical clusters
#' 
#' @description Preserves only one cluster instance per unique element combination
#' 
#' @param c_list Named list with lists of clusters in which to look for duplicates. Assumes the cluster within each named item are all the same size and that cluster elements are all sorted the same ( i.e. all alphabetically sorted ).
#' @param cluster_id_width Width of the cluster id's
#'  
#' @return Returns original list of clusters with duplicates stored in cluster_ids attribute of each unique cluster. Names of list for each size cluster are Cluster_{cluster-#-this_size}x{#-instances-this-cluster} ( i.e. Cluster_001x24 )
#' 
#' @export
#'
collapse_identical_clusters <- function( 
  c_list, 
  cluster_id_width=NA 
){
  cluster_id_width %<>% max(length(strsplit(names(c_list)[length(c_list)], "_")[[1]][[2]]), na.rm=T)
  for (c_index in 1:length(c_list)){
    clust = 1 #which cluster of this size are we looking at?
    clusters_left <- length(c_list[[c_index]]) # how many clusters are there left of this size to compare?
    total_clusters<-clusters_left # how many clusters of this size were there to begin with?
    unique_clusters<-length(unique(c_list[[c_index]])) # how many unique clusters are there?
    while(clust <= clusters_left){ #while there are still clusters to look at
      cluster_to_find <- c_list[[c_index]][[clust]] #what elements are in this cluster?
      #find duplicate clusters
      duplicate_clusters <- which(sapply(c_list[[c_index]], function(x){
        all(x == cluster_to_find)
      }))
      # if any duplicates were found, remove them and set an attribute with list of removed cluster id's in the instance kept
      attr(c_list[[c_index]][[clust]], "cluster_ids") <- names(c_list[[c_index]])[duplicate_clusters]
      c_list[[c_index]][duplicate_clusters[-1]] <- NULL #remove all duplicate clusters from the list ( keeping only the current instance )
      clusters_left <- length(c_list[[c_index]]) # update how many clusters are left since some have probably been removed as duplicates
      clust <- clust+1 #move on to the next cluster of this size
    }
    #sort clusters of this length by frequency of occurrence ( i.e. how many of each unique cluster were found )
    c_list[[c_index]] <- c_list[[c_index]][order( mapply( function(x) length(attr(x,"cluster_ids")), c_list[[c_index]] ), decreasing = T )]
    #preserve total cluster count and unique cluster count for clusters of this size as attributes
    attr(c_list[[c_index]], "total_clusters")<-total_clusters
    attr(c_list[[c_index]], "unique_clusters")<-unique_clusters
    #update list names to Cluster_{#-in-order}x{#-of_duplicates}
    names(c_list[[c_index]]) <- paste(paste("Cluster", stringr::str_pad(1:length(c_list[[c_index]]),width = cluster_id_width, pad = "0"), sep="_"), 
                                      mapply(function(x) length(attr(x,"cluster_ids")),c_list[[c_index]]), sep="x")
  }
  return(c_list)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# print_cluster_counts
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Display cluster distribution within groups by cluster size
#' 
#' @description Prints out elements of unique clusters of each size requested with stats on frequency of occurrence
#' 
#' @param unique_clusters_by_size List of cluster sizes, each size being a list of unique clusters as output by this packages' collapse_identical_clusters method
#' @param cluster_lengths Integer range of cluster sizes to include in output
#'  
#' @return Returns list of requested cluster sizes with elements of unique clusters and stats on frequency
#' 
#' @export
#'
print_cluster_counts <- function( 
  unique_clusters_by_size, 
  cluster_lengths=1:max(mapply(function(x) length(x[[1]]), unique_clusters_by_size)) 
){
  r<-lapply(cluster_lengths, 
            function(x){
              cluster_group <- unique_clusters_by_size[[paste0("clusters_of_",x)]]
              total_clusters <- attr(cluster_group,"total_clusters")
              paste(lapply(names(cluster_group),
                           function(y){
                             cluster_count <- as.integer(strsplit(y,split = "x")[[1]][2])
                             paste(cluster_count,"of",total_clusters, "(",as.integer(cluster_count/total_clusters * 100),"%)")
                           }), 
                    lapply(unique_clusters_by_size[[paste0("clusters_of_",x)]], paste, collapse=", "),
                    sep=" - ")
            })
  names(r)<-paste0("clusters_of_", cluster_lengths)
  r
}