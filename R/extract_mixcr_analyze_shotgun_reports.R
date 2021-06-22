
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# extract_mixcr_analyze_shotgun_reports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Pull the numbers out of the reports from 'mixcr analyze shotgun'
#' 
#' @param input_paths Vector of paths to the mixcr report files
#' 
#' @return Returns data.table with samples by row and mixcr stats by column.
#' 
#' @export
extract_mixcr_analyze_shotgun_reports = function(
  input_paths
){
  for (input_path in input_paths){
    inner_list = list()
    mixcr_report = readr::read_file(input_path)
    mixcr_sections = strsplit(mixcr_report, split="======================================")[[1]]
    run_name = stringr::str_replace(basename(input_path), '.report', '')
    inner_list["Run_ID"] = run_name
    
    # 5 sections: align, ap1, ap2, extend, and assemble
    # these can be distinguished by the output file
    #   HugoLo_IPRES_2016-Pt01-ar-279.vdjca
    #   .rescued_0.vdjca
    #   .rescued_1.vdjca
    #   .extended.vdjca
    #   .clns
    for (mixcr_section in mixcr_sections) {
      section_lines = strsplit(mixcr_section, split="\n")[[1]]
      if (length(section_lines) > 1) {
        output_line = grep("Output file(s):", section_lines, value=T, fixed = T)
        if (grepl(paste0(run_name,".vdjca"), output_line, fixed = T)){
          my_prefix = "Align_"
        } else if (grepl("rescued_0.vdjca$", output_line)) {
          my_prefix = "AP1_"
        } else if (grepl("rescued_1.vdjca$", output_line)) {
          my_prefix = "AP2_"
        } else if (grepl("extended.vdjca$", output_line)) {
          my_prefix = "Extend_"
        } else if (grepl("clns$", output_line)) {
          my_prefix = "Assemble_"
        }
        #message(paste0(my_prefix, "======================================="))
        
        for(section_line in section_lines){
          line_parts = strsplit(section_line, split=": ", fixed = T)[[1]]
          line_name = line_parts[1]
          line_data = strsplit(line_parts[2], split=" (", fixed = T)[[1]][1]
          
          replacement_name = switch(  
            line_name,
            # Align
            "Total sequencing reads" = "Total_Reads",  
            "Successfully aligned reads" = "Aligned Reads",
            "Paired-end alignment conflicts eliminated" = "PE_Conflicts_Eliminated",
            "Alignment failed, no hits (not TCR/IG?)" = "Failed_No_Hits",
            "Alignment failed because of absence of CDR3 parts" = "Failed_No_CDR3",
            "Alignment failed because of low total score" = "Failed_Low_Score",
            "Overlapped" = "Overlapped",
            "Overlapped and aligned" = "Overlapped_Aligned",
            "Alignment-aided overlaps" = "Aided_Overlaps",
            "Overlapped and not aligned" = "Overlapped_Not_Aligned",
            "V gene chimeras" = "V_Chimeras",
            "J gene chimeras" = "J_Chimeras",
            " chains" = "Unknown_Chains",
            "TRA chains" = "TRA_Chains",
            "TRB chains" = "TRB_Chains",
            "TRD chains" = "TRD_Chains",
            "TRG chains" = "TRG_Chains",
            "IGH chains" = "IGH_Chains",
            "TRA,TRD chains" = "TRAD_Chains",
            "TRA chains" = "TRA_Chains",
            "IGK chains" = "IGK_Chains",
            "IGK chains" = "IGK_Chains",
            "IGL chains" = "IGL_Chains",
            "Realigned with forced non-floating bound" = "Realigned",
            "Realigned with forced non-floating right bound in left read" = "Realigned_Right_Bound",
            "Realigned with forced non-floating left bound in right read" = "Realigned_Left_Bound",
            
            # APs
            "Total alignments analysed" = "Total",
            "Number of output alignments" = "Output",
            "Alignments already with CDR3 (no overlapping is performed)" = "No_Overlapping_Needed",
            "Successfully overlapped alignments" = "Overlapping_Needed",
            "Left parts with too small N-region (failed to extract k-mer)" = "Failed_KMer",
            "Extracted k-mer diversity" = "Extracted_KMer",
            "Dropped due to wildcard in k-mer" = "Dropped_Wildcard",
            "Dropped due to too short NRegion parts in overlap" = "Dropped_Short",
            "Dropped overlaps with empty N region due to no complete NDN coverage" = "Dropped_Empty_N",
            "Number of left-side alignments" = "Left_Alignments",
            "Number of right-side alignments" = "Right_Alignments",
            "Complex overlaps" = "Complex_Overlaps",
            "Over-overlaps" = "Over_Overlaps",
            "Partial alignments written to output" = "Output_Partial_Alignments",
            
            # Extension
            "Extended alignments count" = "Alignments",
            "V extensions total" = "V_Total",
            "V extensions with merged targets" = "V_With_Merged",
            "J extensions total" = "J_total",
            "J extensions with merged targets" = "J_With_Merged",
            "V+J extensions" = "VJ_Extensions",
            "Mean V extension length" = "V_Extension_Length",
            "Mean J extension length" = "J_Extension_Length",
            
            # Assemble
            "Final clonotype count" = "Clonotypes",
            "Average number of reads per clonotype" = "Reads_per_Clonotype",
            "Reads used in clonotypes, percent of total" = "Reads_Used",
            "Reads used in clonotypes before clustering, percent of total" = "Reads_Used_Pre_Clustering",
            "Number of reads used as a core, percent of used" = "Core_Reads",
            "Mapped low quality reads, percent of used" = "Mapped_Low_Quality_Reads",
            "Reads clustered in PCR error correction, percent of used" = "Reads_Clustered_PCR_Error",
            "Reads pre-clustered due to the similar VJC-lists, percent of used" = "Reads_Preclustered",
            "Reads dropped due to the lack of a clone sequence, percent of total" = "Dropped_Raads_No_Clone",
            "Reads dropped due to low quality, percent of total" = "Dropped_Raads_Low_Quality",
            "Reads dropped due to failed mapping, percent of total" = "Dropped_Reads_Failed_Map",
            "Reads dropped with low quality clones, percent of total" = "Dropped_Raads_Low_Quality_Clone",
            "Clonotypes eliminated by PCR error correction" = "Dropped_Clonotypes_PCR_Error",
            "Clonotypes dropped as low quality" = "Dropped_Clonotypes_Quality",
            "Clonotypes pre-clustered due to the similar VJC-lists" = "",
            ""
          ) 
          
          if(replacement_name == ""){
            #message(paste0("Could not find name for: ", line_name))
          } else {
            line_name = paste0(my_prefix, replacement_name)
            inner_list[line_name] = as.numeric(line_data)
          }
        }
      }
    }
    my_dt = data.table::rbindlist(list(my_dt,inner_list), use.names = T, fill = T)
  }
  return(my_dt)
}