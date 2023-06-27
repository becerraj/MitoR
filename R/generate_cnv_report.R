#' @import openxlsx
#' @title Generation of report
#' @description Generates the excel with the results of the patient's analysis.
#' You will find each CNV detected with it's information such as size, start, end, 
#' association with Franklin and Varsome, etc.
#' @param CNV_calls dataframe generated using AnalyzePatient function.  
#' @param graph ggplot that represents the CNVs.
#' @param path_dir path of the directory of the patient that has been analized. It must contain R1 and R2.
#' @param transit transition probability: Transition probability of the hidden Markov Chain 
#' from the normal copy number state to either a deletion or a duplication. Is the value used
#' in the generation of the CNV_calls indicated in the first parameter. 
#' @return path where the report has been saved.
#' @export

generate_cnv_report <- function(CNV_calls, graph, path_dir, transit) {
  
  id <- basename(path_dir)
  
  wb <- openxlsx::createWorkbook()
  
  COLNAMES_STYLE <- createStyle(
    fontSize = 12,
    textDecoration = "bold",
    halign = "center", valign = "center", border = "TopBottom",
    borderColour = "black")
  
  
  hyperlink_style <- openxlsx::createStyle(fontColour = "#0000FF")
  
  addWorksheet(wb, paste(id, "CNVAnalisis"))
  addWorksheet(wb, paste(id, "CNVplots"))
  writeData(wb, sheet = 1, x = CNV_calls,
            headerStyle = COLNAMES_STYLE,
            borderStyle = "dashed",
            borders = "columns", borderColour = "black")
  
  
  setColWidths(wb, sheet = 1, cols = 1:4, widths = "auto")
  setColWidths(wb, sheet = 1, cols = 6:ncol(CNV_calls), widths = "auto")
  setColWidths(wb, sheet = 1, cols = 5, widths = 16)
  
  for(i in 1:nrow(CNV_calls)){
    x <- CNV_calls$Franklin[i]
    y <- CNV_calls$Varsome[i]
    names(x) <- c("View on Franklin DB")
    names(y) <- c("View on Varsome DB")
    class(x) <- "hyperlink"
    class(y) <- "hyperlink"
    
    indiceF <- which(colnames(CNV_calls) == "Franklin")
    indiceV <- which(colnames(CNV_calls) == "Varsome")
    
    writeData(wb, sheet = 1, x = x, startRow= i+1, startCol = indiceF)
    writeData(wb, sheet = 1, x = y, startRow= i+1, startCol = indiceV)
    
  }
  
  #Colour cells if the size of the variation is >500
  
  for (i in 1:nrow(CNV_calls)){
    if (CNV_calls$size[i] > 500){
      addStyle(wb, 1, style = openxlsx::createStyle(fgFill = "thistle"), rows = i+1, cols = 1:ncol(CNV_calls))
    }
  }
  
  #Add Plot
  print(graph)
  insertPlot(wb, 2, xy = c(1, 1), width = 25, height = 18, fileType = "png", units = "cm")
  
  openxlsx::saveWorkbook(wb, paste(path_dir, "/MitoR/", id, "_CNVreport_", transit, ".xlsx", sep = ""),  overwrite = TRUE)
  browseURL(paste(path_dir, "/MitoR/", id, "_CNVreport_", transit, ".xlsx", sep = ""))
  
  return(paste(dir_patient, "/MitoR/", id, "_report_", transit, ".xlsx", sep = ""))
  
}
