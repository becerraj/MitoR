#' @title Open
#' @description MitoR-.
#' @import openxlsx

open_MitoR_DataBase <- function(){
  #mitor_files <- sprintf("%s/mitor", Sys.getenv('R_LIBS_USER'))
  mitor_files <- "/home/juan"
  if (!file.exists(sprintf("%s/MitoRSoftware/RDS_DB.rds", mitor_files))){
    stop("There are not any patients on the MitoR DataBase
         In order to see see your DataBase, you must have analyzed at least one patient first.")
  }
  RDS_DB <- readRDS(sprintf("%s/MitoRSoftware/RDS_DB.rds", mitor_files))

  MitoR_Stats <- RDS_DB$Freq_MitoR[[1]]

  Filter_Params <- sapply(RDS_DB[2:length(RDS_DB)], function(x) x[[2]])
  Analysis_Date <- sapply(RDS_DB[2:length(RDS_DB)], function(x) x[[3]])
  Info <- data.frame(Patients = names(RDS_DB[2:length(RDS_DB)]), Filter_Params = Filter_Params, Analysis_Date = Analysis_Date)

  ### GENERACION DE XLSX ###
  #  Crea el workbook
  #wb <- openxlsx::createWorkbook()
  wb <- createWorkbook()
  
  # Estilos de las fuentes
  COLNAMES_STYLE <- openxlsx::createStyle(
    fontSize = 12,
    textDecoration = "bold",
    halign = "center", valign = "center", border = "TopBottom",
    borderColour = "black")

  hyperlink_style <- openxlsx::createStyle(fontColour = "#0000FF")

  # Sheets dentro del excel
  openxlsx::addWorksheet(wb, "MitoR_Stats")
  openxlsx::addWorksheet(wb, "Info")

  # Agregamos el dataframe y su diseno
  openxlsx::writeData(wb, sheet = 1, x = MitoR_Stats,
            headerStyle = COLNAMES_STYLE,
            borderStyle = "dashed",
            borders = "columns", borderColour = "black")

  openxlsx::writeData(wb, sheet = 2, x = Info,
                      headerStyle = COLNAMES_STYLE,
                      borderStyle = "dashed",
                      borders = "columns", borderColour = "black")


  # Reemplazamos los valores porel hipervinculo
  for (i in 1:nrow(MitoR_Stats)){
    x <- MitoR_Stats$Franklin[i]
    names(x) <- c("View on Franklin DB")
    class(x) <- "hyperlink"
    y <- MitoR_Stats$VarSome[i]
    names(y) <- c("View on VarSome DB")
    class(y) <- "hyperlink"
    z <- MitoR_Stats$dbSNP[i]
    names(z) <- c("View on dbSNP DB")
    class(z) <- "hyperlink"

    # writeDara sirve para escribir datos sobre el wb. En nuestro caso queremos cambiar los hipervinculos aca
    openxlsx::writeData(wb, sheet = 1, x = x, startRow = i+1, startCol = 13)
    openxlsx::writeData(wb, sheet = 1, x = y, startRow = i+1, startCol = 14)
    openxlsx::writeData(wb, sheet = 1, x = z, startRow = i+1, startCol = 9)
  }
  # Crea un archivo temporal en la memoria y escribe el dataframe en él
  archivo_temporal <- tempfile(fileext = ".xlsx")
  openxlsx::saveWorkbook(wb, file = archivo_temporal, overwrite = TRUE)

  # Abre el archivo temporal en la aplicación asociada a los archivos XLSX
  system(paste("open", archivo_temporal))
  # browseURL(archivo_temporal)
  # system(sprintf('start %s', archivo_temporal), wait = FALSE)
}
