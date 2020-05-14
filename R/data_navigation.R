
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' get_features_from_gmt_file
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title get_features_from_gmt_file 
#' 
#' @description 
#' Just grabs the names of all of the gene sigs from a gmt file.
#' 
#' @param gmt_path	A string path to a gmt text file.
#' 
#' @return A vector strings which includes the gene signature names 
#' 
#' @family gene_signatures
#' 
#' @export
get_features_from_gmt_file = function(gmt_path){
  lapply(readLines(gmt_path), function(x){
    my_split = strsplit(x, "\t") %>% unlist
    my_split[1]
  }) %>% unlist()
}



#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' get_gene_set_entrez_human_bgvlab_path
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title get_gene_set_entrez_human_bgvlab_path 
#' 
#' @description 
#' Returns the path to the entrez_ids, human, BGVLab signatures in the package library.
#' 
#' @param none
#' 
#' @return A path to the gmt file.
#' 
#' @family gene_signatures
#' 
#' @export
get_gene_set_entrez_human_bgvlab_path = function(){
  return(system.file(file.path("gene_sets", "entrez_ids", "human", "BGVLab", "BGVLab_entrez_human.gmt.txt"), package = "binfotron"))
  }


#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' get_gene_set_hgnc_human_bgvlab_path
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title get_gene_set_hgnc_human_bgvlab_path 
#' 
#' @description 
#' Returns the path to the symbol, human, BGVLab signatures in the package library.
#' 
#' @param none
#' 
#' @return A path to the gmt file.
#' 
#' @family gene_signatures
#' 
#' @export
get_gene_set_hgnc_human_bgvlab_path = function(){
  return(system.file(file.path("gene_sets", "symbol", "human", "BGVLab", "hgnc_human_bgvlab.gmt.txt"), package = "binfotron"))
}


#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' get_gene_set_entrez_human_kegg_path
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title get_gene_set_entrez_human_kegg_path 
#' 
#' @description 
#' Returns the path to the entrez_id, human, kegg signatures in the package library.
#' 
#' @param none
#' 
#' @return A path to the gmt file.
#' 
#' @family gene_signatures
#' 
#' @export
get_gene_set_entrez_human_kegg_path = function(){
  return(system.file("gene_sets", "entrez_ids", "human", "KEGG", "kegg_hsa_entrez.gmt.txt", package = "binfotron"))
}



#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' get_gene_set_hgnc_human_kegg_path
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title get_gene_set_hgnc_human_kegg_path 
#' 
#' @description 
#' Returns the path to the human genome nomenclature, kegg signatures in the package library.
#' 
#' @param none
#' 
#' @return A path to the gmt file.
#' 
#' @family gene_signatures
#' 
#' @export
get_gene_set_hgnc_human_kegg_path = function(){
  return(system.file("gene_sets", "symbol", "human", "KEGG", "kegg_hsa_hgnc.gmt.txt", package = "binfotron"))
}


#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' get_gene_set_entrez_human_c2_path
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title get_gene_set_entrez_human_c2_path 
#' 
#' @description 
#' Returns the path to the entrez_id, human, c2 signatures in the package library.
#' 
#' @param none
#' 
#' @return A path to the human_c2_v5p1.rdata file.
#' 
#' @family gene_signatures
#' 
#' @export
get_gene_set_entrez_human_c2_path = function(){
  
  return(system.file("gene_sets", "entrez_ids", "human", "MSigDB", "human_c2_v5p1.rdata", package = "binfotron"))
}


#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' get_gene_set_hgnc_human_c2cp_path
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title get_gene_set_hgnc_human_c2cp_path 
#' 
#' @description 
#' Returns the path to the entrez_id, human, c2 signatures in the package library.
#' 
#' @param none
#' 
#' @return A path to the human_c2_v5p1.rdata file.
#' 
#' @family gene_signatures
#' 
#' @export
get_gene_set_hgnc_human_c2cp_path = function(){
  return(system.file("gene_sets", "symbol", "human", "c2cp", "c2.cp.v6.0.symbols.gmt.txt", package = "binfotron"))
}



#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' get_biomart_hsa_ucsc_path
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title get_biomart_hsa_ucsc_path 
#' 
#' @description 
#' Returns the path a saved mart
#' 
#' @param none
#' 
#' @return A path to the rds file.
#' 
#' @family mart
#' 
#' @export
get_biomart_hsa_ucsc_path = function(){
  return(system.file("biomart", "hsa_ensembl_ucsc", "hsa_ensembl_ucsc.rds", package = "binfotron"))
}


#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' get_corrected_entropy_rdata_2_4096_path
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title get_corrected_entropy_rdata_2_4096_path 
#' 
#' @description 
#' Returns the path to the corrected entropy model for abundances in the range of 2-4096.
#' 
#' @return A path to the RData file
#' 
#' @family diversity
#' 
#' @export
get_corrected_entropy_rdata_2_4096_path = function(){
  return(system.file("models", "subsampled_stig_2-4096_model", 
                               "RNASeq_subsampled_random_abundance_TRB_optimal_model.RData", 
                               package = "binfotron"))
}


#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' get_corrected_entropy_rdata_4097_65536_path
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title get_corrected_entropy_rdata_2_65536_path 
#' 
#' @description 
#' Returns the path to the corrected entropy model for abundances in the range of 4097-65536.
#' 
#' @return A path to the RData file
#' 
#' @family diversity
#' 
#' @export
get_corrected_entropy_rdata_4097_65536_path = function(){
  return(system.file("models", "subsampled_stig_4097-65536_model", 
                               "RNASeq_subsampled_random_abundance_TRB_optimal_model.RData", 
                               package = "binfotron"))
}




#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' get_corrected_entropy_rdata_ab8_1024_ent1_8_path
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title get_corrected_entropy_rdata_ab8_1024_ent1_8_path 
#' 
#' @description 
#' Returns the path to the corrected entropy model for abundances in the range of 1025-65536.
#' 
#' @return A path to the RData file
#' 
#' @family diversity
#' 
#' @export
get_corrected_entropy_rdata_ab8_1024_ent1_8_path = function(){
  return(system.file("models", "tcga_ab8-1024_ent1-8", 
                     "tcga_distribution_optimal_model.RData", 
                     package = "binfotron"))
}

#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' get_corrected_entropy_rdata_ab1025_65K_ent1-8_path
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title get_corrected_entropy_rdata_ab1025_65K_ent1 
#' 
#' @description 
#' Returns the path to the corrected entropy model for abundances in the range of 1025-65536.
#' 
#' @return A path to the RData file
#' 
#' @family diversity
#' 
#' @export
get_corrected_entropy_rdata_ab1025_65K_ent1_8_path = function(){
  return(system.file("models", "tcga_ab1025-65K_ent1-8", 
                     "tcga_distribution_optimal_model.RData", 
                     package = "binfotron"))
}



#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' get_gene_set_entrez_human_cta_path
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title get_gene_set_entrez_human_cta_path 
#' 
#' @description 
#' Returns the path to the entrez_ids, human, CTA signatures in the package library.
#' 
#' @param none
#' 
#' @return A path to the gmt file.
#' 
#' @family gene_signatures
#' 
#' @export
get_gene_set_entrez_human_cta_path = function(){
  return(system.file(file.path("gene_sets", "entrez_ids", "human", "CTA", "CTA_entrez.gmt.txt"), package = "binfotron"))
}

#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' get_gene_set_hgnc_human_cta_path
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title get_gene_set_hgnc_human_cta_path 
#' 
#' @description 
#' Returns the path to the hgnc, human, CTA signatures in the package library.
#' 
#' @param none
#' 
#' @return A path to the gmt file.
#' 
#' @family gene_signatures
#' 
#' @export
get_gene_set_hgnc_human_cta_path = function(){
  return(system.file(file.path("gene_sets", "symbol", "human", "cta", "CTA_hgnc.gmt.txt"), package = "binfotron"))
}
