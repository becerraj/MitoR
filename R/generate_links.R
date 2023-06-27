#' @title Generation of Links
#' @description Generates links to databases for each variation
#' @param CNV_calls dataframe created with ExomeDepth algorithm.
#' @return CNV_calls with 2 new columns "Franklin" and "Varsome" that contain
#' the links that can be pasted in the browser to be redirected to any of 
#' these databases and get all the information about t¡an especific variation
#' within the CNVs found in this case.

generate_links <- function(CNV_calls) {
  
  for (i in 1:nrow(CNV_calls)) {
    mutF <- CNV_calls$id[i]
    mutF <- gsub(':', '-', mutF)
    mutV <- gsub("-", "%3A", mutF)
    addressF <- paste("https://franklin.genoox.com/clinical-db/variant/sv/", mutF, sep = "")
    addressV <- paste("https://varsome.com/cnv/hg19/", mutV, sep = "")
    
    if (CNV_calls$type[1] == "deletion") { 
      addressF <- paste(addressF, "-DEL", sep = "")
      addressV <- paste(addressV, "%3ADEL?", sep = "")
    } else {
      addressF <- paste(addressF, "-DUP", sep = "")
      addressV <- paste(addressV, "%3ADUP?", sep = "")
    }
    
    CNV_calls$Franklin[i] <- addressF
    CNV_calls$Varsome[i] <- addressV
    
  }
  
  CNV_calls <- CNV_calls[, c(1,2,3,4,5,6,11,12,7,8,9,10)]
  return(CNV_calls)
}