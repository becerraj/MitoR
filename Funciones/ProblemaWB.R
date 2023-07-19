VCF_file <- "/home/daniela/Paciente/4456999/MitoR/4456999_filtered_MitoR.vcf"


wb_list <- generateXLSX(VCF_file)
wb_list[[3]] <- xlsx_plot_mut_sheet

wb_complete <- openxlsx::createWorkbook()

for (wb_actual in wb_list) {
  num_sheets <- length(wb_actual$worksheets)
  cs <- length(wb_complete$worksheets) +1
  for (s in 1:num_sheets) {
    print(s)
    datos <- openxlsx::read.xlsx(wb_actual, sheet = s)
    addWorksheet(wb_complete, sheetName = cs)
    writeData(wb_complete, sheet = cs, x = datos)
    cs <- cs + 1
  }
}

openxlsx::saveWorkbook(wb_complete, "/home/daniela/prueba_wb.xlsx",  overwrite = TRUE)
browseURL("/home/daniela/prueba_wb.xlsx")

openxlsx::saveWorkbook(xlsx_plot_mut_sheet, "/home/daniela/graficos.xlsx",  overwrite = TRUE)
browseURL("/home/daniela/graficos.xlsx")




#----------------------
wb1 <- openxlsx::createWorkbook()
wb2 <- openxlsx::createWorkbook()


COLNAMES_STYLE <- createStyle(
  fontSize = 12,
  textDecoration = "bold",
  halign = "center", valign = "center", border = "TopBottom",
  borderColour = "black")


addWorksheet(wb1, "hoja1")
addWorksheet(wb1, "hoja2")
writeData(wb1, sheet = 1, x = CNV_calls,
          headerStyle = COLNAMES_STYLE,
          borderStyle = "dashed",
          borders = "columns", borderColour = "black")
writeData(wb1, sheet = 2, x = CNV_calls,
          headerStyle = COLNAMES_STYLE,
          borderStyle = "dashed",
          borders = "columns", borderColour = "black")

addWorksheet(wb2, "hoja1")
writeData(wb2, sheet = 1, x = PatientsDB,
          headerStyle = COLNAMES_STYLE,
          borderStyle = "dashed",
          borders = "columns", borderColour = "black")



wb_list<-list()
wb_list[[1]]<-wb1
wb_list[[2]]<-wb2

wb_complete <- openxlsx::createWorkbook()

for (wb_actual in wb_list) {
  num_sheets <- length(wb_actual$worksheets)
  cs <- length(wb_complete$worksheets) +1
  for (s in 1:num_sheets) {
    print(s)
    datos <- openxlsx::read.xlsx(wb_actual, sheet = s)
    addWorksheet(wb_complete, sheetName = cs)
    writeData(wb_complete, sheet = cs, x = datos)
    cs <- cs + 1
  }
}


openxlsx::saveWorkbook(wb_complete, "/home/daniela/prueba_wb.xlsx",  overwrite = TRUE)
browseURL("/home/daniela/prueba_wb.xlsx")

