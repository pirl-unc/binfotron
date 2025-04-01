#' encoding_decoding
#' there are a lot of characters that aren't friendly for data.frame,/data.table
#'   operations as well as modeling and other functions.  Functions here will 
#'   encode r decode those values so they don't need so much manual intervention

#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' encode_special_char
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' @title encodes non alpha numeric or underscore characters into hex.
#' 
#' @description
#' Non alpha numeric text (eg [a-zA-Z0-9_]) is converted to its unicode hex value
#' _hx[0-9A-F]{2}_
#' 
#' @param char_vector character vector that will be encoded
#' 
#' @export
encode_special_char = function(char_vector) {
	sapply(char_vector, function(name) {
		encoded_name = ""
		name = gsub(" ", "_", name)
		for (char in strsplit(name, '')[[1]]) {
			if (grepl("[a-zA-Z0-9_]", char)) {
				encoded_name = paste0(encoded_name, char)
			} else {
				encoded_name = paste0(encoded_name, "_hx", toupper(as.character(sprintf("%x", as.integer(charToRaw(char))))), "_")
			}
		}
		return(encoded_name)
	})
}


#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' decode_special_char
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' @title Decodes unicode hex value
#' 
#' @description
#' Unicode hex value _hx[0-9A-F]{2}_ is converted into non alpha numeric text
#' 
#' @param char_vector character vector that will be decoded
#' 
#' @export
decode_special_char = function(char_vector) {
	sapply(char_vector, function(name) {
		while (grepl("_hx[0-9A-F]{2}_", name)) {
			my_match = regmatches(name, regexpr("_hx[0-9A-F]{2}_", name))
			if(length(my_match) > 0){
				my_match = my_match[[1]]
			  hex_code = gsub("_hx|_", "", my_match)
			  char = rawToChar(as.raw(strtoi(hex_code, base = 16L)))
			  name = sub(my_match, char, name, fixed = TRUE)
			}
		}
	  name = gsub("_", " ", name)
		return(name)
	})
}



#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' decode_clms
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' @title Decodes data.frame or data.table using decode_special_char
#' 
#' @description
#' Unicode hex value _hx[0-9A-F]{2}_ columns are converted into non alpha numeric text
#' 
#' @param my_df Input date.frame or data.table
#' @param skip_clms character vector of columns to skip
#' @param replace_existing_clm Boolean to indicate whether encoded columns should be overwritten or written as <clm_name>_Decoded
#' 
#' @export
decode_clms = function( my_df, skip_clms = c(), replace_existing_clm = FALSE){
	for (clm_name in names(my_df)) {
		if ( ! clm_name %in% skip_clms ){
			clm_values = my_df[[clm_name]]
			if (any(grepl("_hx[0-9A-F]{2}_", clm_values))) {
				if (replace_existing_clm){
					decoded_clm_name = clm_name
				} else {
					decoded_clm_name = paste0(clm_name, "_Decoded")
				}
				my_df[[decoded_clm_name]] = decode_special_char(clm_values)
			}
		}
	}
	return(my_df)
}