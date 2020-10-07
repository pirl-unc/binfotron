
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' create_gene_lookup
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' @title Create gene lookup 
#' 
#' @description
#' Uses \code{\link[AnnotationDbi]{select}} to produce a lookup array for converting gene names
#' If using the default value for \code{db_object}(\code{org.Hs.eg.db}) the values for the inputs 
#' and outputs include any listed in \url{https://www.bioconductor.org/help/course-materials/2014/useR2014/Integration.html},
#' although the most common will likely be: ENTREZID, ACCNUM, ALIAS, UNIGENE, ENSEMBL, ENSEMBLPROT, 
#' ENSEMBLTRANS, GENENAME, UNIPROT, OMIM, UCSCKG, SYMBOL
#' 
#' @param sep String to indicate how gene names should be collapsed.  a value of 
#' \code{NULL} will return just the first item of the returned results
#' @param input_values Character vecter of names to lookup
#' @param input_type String of the input values. See description for common options.
#' @param output_type String of the ouput values. See description for common options.
#' @param db_object AnnotationDb object.  If left null it will be set to \code{org.Hs.eg.db}.
#' 
#' @return Named array that can be used to lookup gene values
#' 
#' @export
create_gene_lookup = function(
  input_values,
  input_type,
  output_type,
  sep = "\t",
  db_object = NULL
){
  # TODOs add cross species lookups
  
  library(annotate)
  
  if(is.null(db_object)){
    library(org.Hs.eg.db)
    db_object = org.Hs.eg.db
  }
  return_values = AnnotationDbi::select(org.Hs.eg.db, keys=input_values, columns=output_type, keytype=input_type)
  fwd_counts = tapply(return_values[[2]], return_values[[1]], function(x){length(x)})
  rev_counts = tapply(return_values[[1]], return_values[[2]], function(x){length(x)})
  # report:
  #  1 index to mean X outputs
  fwd_mean = mean(fwd_counts, na.rm = TRUE)
  rev_mean = mean(rev_counts, na.rm = TRUE)
  
  cat(paste0("Each input mapped to average of ", fwd_mean," output(s)", "\n"))
  cat(paste0("Each output mapped to average of ", rev_mean," input(s)", "\n"))
  
  if(is.null(sep)){
    output_values = tapply(return_values[[2]], return_values[[1]], function(x){
      x = x[!is.na(x)]
      return(x[1])
    })
  } else {
    output_values = tapply(return_values[[2]], return_values[[1]], function(x){paste0(x[!is.na(x)], collapse = sep)})
  }
  
  return(output_values)
}



#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' convert_gmt_file
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' @title Create gene lookup 
#' 
#' @description
#' Uses \code{\link{create_gene_lookup}} to convert one type of gmt file to another.
#' If using the default value for \code{db_object}(\code{org.Hs.eg.db}) the values for the inputs 
#' and outputs include any listed in \url{https://www.bioconductor.org/help/course-materials/2014/useR2014/Integration.html},
#' although the most common will likely be: ENTREZID, ACCNUM, ALIAS, UNIGENE, ENSEMBL, ENSEMBLPROT, 
#' ENSEMBLTRANS, GENENAME, UNIPROT, OMIM, UCSCKG, SYMBOL
#' 
#' @param input_path Path to the gmt file you would like ot convert.
#' @param input_type Typically ENTREZID or SYMBOL. See decription.
#' @param output_type Typically SYMBOL or ENTREZID. See decription. 
#' @param output_path Full path to output gmt file.
#' @param db_object AnnotationDb object. If left null it will be set to  \code{org.Hs.eg.db}. Passed to 
#' \code{\link{create_gene_lookup}}.
#' @param na_string String to fill in the na positions of a gmt file. Typically "NA" or "".
#' @param input_readme_path Optional path to readme to import into this gmt's readme.
#'
#' @return None
#' 
#' @export
convert_gmt_file = function(
  input_path,
  input_type,
  output_type,
  output_path,
  db_object = NULL,
  na_string = "NA",
  input_readme_path = NULL
){
  if(is.null(db_object)){
    library(org.Hs.eg.db)
    db_object = org.Hs.eg.db
  }
  
  gmt_lines = readLines(input_path)
  all_genes = c()
  combined_list = list()
  max_sig_length = 0
  for(sig_index in 1:length(gmt_lines)){
    gmt_line = gmt_lines[sig_index]
    gmt_split = strsplit(gmt_line, split = "\t")[[1]]
    sig_list = list(ref = gmt_split[2], genes = gmt_split[3:length(gmt_split)])
    sig_list$genes = unique(sig_list$genes[sig_list$genes != "NA" & sig_list$genes != ""])
    if(length(sig_list$gene) > max_sig_length) max_sig_length = length(sig_list$gene)
    all_genes = c(all_genes, sig_list$genes)
    combined_list[[gmt_split[1]]] = sig_list
  }
  all_genes = unique(all_genes)
  
  gene_lu = create_gene_lookup(
    input_values= all_genes,
    input_type = input_type,
    output_type = output_type,
    sep = "\t",
    db_object = org.Hs.eg.db
  )
  
  for(sig_index in 1:length(combined_list)){
    my_genes = gene_lu[combined_list[[sig_index]]$genes]
    na_needed = max_sig_length - length(my_genes)
    if(na_needed > 0) my_genes = c(my_genes, rep(na_string, na_needed))
    combined_list[[sig_index]]$genes = paste0(my_genes, collapse = "\t")
  }
  
  if(is.null(input_readme_path)){
    input_readme_path = list.files(dirname(input_path), pattern = "readme", ignore.case = TRUE, full.names = TRUE)[1]
  }
  
  readme_notes = paste0("Started with gmt file: ", input_path)
  readme_notes = c(readme_notes, "")
  if(!is.na(input_readme_path)){
    readme_notes = c(readme_notes, housekeeping::import_annotation(input_readme_path))
  }
  readme_notes = c(readme_notes, 
                      paste0("On ", format(Sys.time(), "%a %b %d %X %Y"), ", converted the gmt file from ", 
                             input_type, " to ", output_type, 
                             " using binfotron::convert_gmt_file v", packageVersion("binfotron")))
  
  output_dir = dirname(output_path)
  dir.create(output_dir, showWarnings = F)
  writeLines(readme_notes, file.path(output_dir, "readme.txt"))
  
  if( file.exists(output_path)) unlink(output_path)
  for(sig_index in 1:length(combined_list)){
    sig_name = names(combined_list)[sig_index]
    sig_ref = combined_list[[sig_index]]$ref
    sig_genes = combined_list[[sig_index]]$genes
    write(paste0(sig_name, "\t", sig_ref, "\t", sig_genes),output_path, append = TRUE)
  }
}
