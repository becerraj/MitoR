#' @title View Data Bases
#' @description Download the databases in excel file.
#' @export
#' @import openxlsx

view_dbs <- function() {

  mutDB <- readRDS(paste(Sys.getenv('R_LIBS_USER'), "/mitorDB/RDS_DB.RDS", sep=""))
  mut_DB <- mutDB[[1]]

  cnv_DBs <- readRDS(paste(Sys.getenv('R_LIBS_USER'), "/mitorDB/cnvDB.RDS", sep=""))
  ControlDB <- cnv_DBs[[1]]
  PatientsDB <- cnv_DBs[[2]]
  GenesDB <- cnv_DBs[[3]]

  wb_dbs <- openxlsx::createWorkbook()

  COLNAMES_STYLE <- createStyle(
    fontSize = 12,
    textDecoration = "bold",
    halign = "center", valign = "center", border = "TopBottom",
    borderColour = "black")

  hyperlink_style <- openxlsx::createStyle(fontColour = "#0000FF")

  addWorksheet(wb_dbs, "MutationsDB")
  addWorksheet(wb_dbs, "GenesDB")
  addWorksheet(wb_dbs, "PatientsDB")
  addWorksheet(wb_dbs, "ControlDB")

  writeData(wb_dbs, sheet = 1, x = mut_DB,
            headerStyle = COLNAMES_STYLE,
            borderStyle = "dashed",
            borders = "columns", borderColour = "black")
  writeData(wb_dbs, sheet = 2, x = GenesDB,
            headerStyle = COLNAMES_STYLE,
            borderStyle = "dashed",
            borders = "columns", borderColour = "black")
  writeData(wb_dbs, sheet = 3, x = PatientsDB,
            headerStyle = COLNAMES_STYLE,
            borderStyle = "dashed",
            borders = "columns", borderColour = "black")
  writeData(wb_dbs, sheet = 4, x = ControlDB,
            headerStyle = COLNAMES_STYLE,
            borderStyle = "dashed",
            borders = "columns", borderColour = "black")


  saveWorkbook(wb_dbs, paste(Sys.getenv('R_LIBS_USER'), "/mitorDB/DataBases.xlsx", sep=""), overwrite = TRUE)
  browseURL(paste(Sys.getenv('R_LIBS_USER'), "/mitorDB/DataBases.xlsx", sep=""))

}
