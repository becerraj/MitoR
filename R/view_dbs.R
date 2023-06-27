#' @title View Data Bases
#' @description Download the databases in excel file.
#' @export
#' @import openxlsx

view_dbs <- function() {
  
  wb <- openxlsx::createWorkbook()
  
  COLNAMES_STYLE <- createStyle(
    fontSize = 12,
    textDecoration = "bold",
    halign = "center", valign = "center", border = "TopBottom",
    borderColour = "black")
  
  hyperlink_style <- openxlsx::createStyle(fontColour = "#0000FF")
  
  addWorksheet(wb, "GenesDB")
  addWorksheet(wb, "PatientsDB")
  addWorksheet(wb, "ControlDB")
  
  initial <- initialize_db()
  ControlDB <- initial[[1]]
  PatientsDB <- initial[[2]]
  GenesDB <- initial[[3]]
  
  writeData(wb, sheet = 1, x = GenesDB,
            headerStyle = COLNAMES_STYLE,
            borderStyle = "dashed",
            borders = "columns", borderColour = "black")
  writeData(wb, sheet = 2, x = PatientsDB,
            headerStyle = COLNAMES_STYLE,
            borderStyle = "dashed",
            borders = "columns", borderColour = "black")
  writeData(wb, sheet = 3, x = ControlDB,
            headerStyle = COLNAMES_STYLE,
            borderStyle = "dashed",
            borders = "columns", borderColour = "black")
  

  saveWorkbook(wb, paste(Sys.getenv('R_LIBS_USER'), "/mitorDB/DataBases.xlsx", sep=""), overwrite = TRUE)
  browseURL("~/DataBases.xlsx")
  
}