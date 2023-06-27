#' @title Update the MitoR-DataBase stats taken from HMTVAR
#' @description MitoR-DataBase takes information for each mutation from the HMTVAR DataBase. This DataBase is constantly updated.
#' We recommend to update the information of MitoR-DataBase at least once a month. With this function, the MitoR DataBase will be updated with the last
#' available information about each mutation from HMTVAR.
#' @param AddToPatients In case you want to update the information of the individual XLSX files from the patients as well, keep AddToPatients = TRUE.
#' Please be aware that this will take a while if you have many patients to update.
#' In case you want to update only one patient or just some of them, we recommend you to create the XLSX file from the VCF again. Check generateXLSX() for more information.
#' @return The DataBase XLSX file will automatically be updated. The same will happen with the patients XLSX files in case of AddToPatients = TRUE
#' @examples
#' updateMitoRStats_HMTVAR()
#' updateMitoRStats_HMTVAR(addToPatients = TRUE)

updateMitoRStats_HMTVAR <- function(AddToPatients = FALSE){
  # Actualiza los datos de la DB de HMTVAR
  mitor_files <- sprintf("%s/mitor", Sys.getenv('R_LIBS_USER'))
  RDS_DB <- readRDS(sprintf("%s/MitoRSoftware/RDS_DB.rds", mitor_files))

  Freq_MitoR <- RDS_DB$Freq_MitoR[[1]]

  # La funcion buscarHMTVAR hace las busquedas en base al DataFrame que le damos
  # Actualizamos el DataFrame del MitoRStats
  updates <- buscarHMTVAR(Freq_MitoR)
  Freq_MitoR$clinVAR <- updates[[1]]
  Freq_MitoR$MitoMap <- updates[[2]]
  Freq_MitoR$dbSNP <- updates[[3]]
  Freq_MitoR$Omim <- updates[[4]]
  Freq_MitoR$Disease <- updates[[5]]
  Freq_MitoR$Freq_HMTVAR <- updates[[6]]

  # Deberia hacer un listado de ultimas actualizaciones de la base de datos y que se agregue en la segunda hoja

  # Si no es la primera actualizacion, borramos la fila que marca la actualizacion anterior
  updateRow <- data.frame(Paciente = "HMTVAR Update", Date = format(Sys.Date(), "%d/%m/%Y"), version_BWA = "-", version_GATK = "-", version_PICARD = "-", Last_Update = "-")
  if ("HMTVAR Update" %in% SoftData$Paciente) {
    SoftData <- SoftData[!SoftData$Paciente == "HMTVAR Update", ]
  }

  # Agregamos la fila de la actualizacion
  SoftData <<- rbind(SoftData, updateRow)
  openxlsx::write.xlsx(mitoRStats, MitoRStats_Path, sheetName = "MitoRStats" , rowNames = FALSE)
  openxlsx::write.xlsx(SoftData, MitoRStats_Path, sheetName = "SoftData", rowNames = FALSE, overwrite = TRUE)

  if (AddToPatients == TRUE) {
    for (i in names(RDS_DB)[2:length(RDS_DB)]) {
      tryCatch(
        expr = {
          updateFreqMutMitoR(i)
        },
        error = function(e) {
          message(sprintf("An error occured while updating the patient number '%s'.
                          Please, try to update it individually by using the updateFreqMutMitoR() function.", i))

          print(e)
        },
        warning = function(w) {
          message(sprintf("An error occured while updating the patient number '%s'.
                          Please, try to update it individually by using the updateFreqMutMitoR() function.", i))
          print(w)
        },
        finally = {
          message(sprintf("Patient number '%s' has been succesfully updated.", i))
        }
      )
    }
  }
  RDS_DB$Freq_MitoR <- Freq_MitoR
  saveRDS(RDS_DB, sprintf("%s/MitoRSoftware/RDS_DB.rds", mitor_files))
}

#' @title Update the MitoR-DataBase frequency
#' @description MitoR-
#' @param numPaciente blabla
#' @import openxlsx

updateFreqMutMitoR <- function(numPaciente){ # Agrega o actualiza las frecuencias de la MitoRStats en los XLSX de los pacientes
  mitor_files <- sprintf("%s/mitor", Sys.getenv('R_LIBS_USER'))
  RDS_DB <- readRDS(sprintf("%s/MitoRSoftware/RDS_DB.rds", mitor_files))

  Freq_MitoR <- RDS_DB$Freq_MitoR[[1]]

  patient_XLSX_file <- RDS_DB[[numPaciente]]["XLSX_Path"]
  # Lee el xlsx del paciente: 1 - mutaciones. 2 - softwares
  # A uno le agregamos la frecuencia de la DB. Al otro le agregamos fecha de update
  paciente <- openxlsx::read.xlsx(patient_XLSX_file, sheet = 1)
  softwares <- openxlsx::read.xlsx(patient_XLSX_file, sheet = 2)

  # Tomamos las mutaciones que tiene el paciente
  mutations <- RDS_DB[[numPaciente]]["Mutations"]
  toAddFreq <- strsplit(unlist(mutations, use.names = FALSE), "/")

  # Guardamos en una variable nueva las frecuencias que hay en la DB sobre las mutaciones del paciente
  freqMut_MitoRStats <- c()
  for (i in (1:length(toAddFreq))){
    freqMut_MitoRStats <- c(freqMut_MitoRStats, Freq_MitoR[(Freq_MitoR$POS == as.integer(toAddFreq[[i]][2])) & (Freq_MitoR$ALT == toAddFreq[[i]][3]), ]$Freq_MitoR)
  }

  if ("Freq_MitoR" %in% colnames(paciente)){ # Si es actualizacion de la frecuencia
    paciente$Freq_MitoR <- freqMut_MitoRStats
    softwares$Last_Upd_FreqMitoR <- format(Sys.Date(), "%d/%m/%Y")
  } else{ # Si es la primera vez que agregamos la frecuencia de la base de datos
    paciente <- cbind(paciente, Freq_MitoR = freqMut_MitoRStats)
    softwares <- cbind(softwares, Last_Upd_FreqMitoR = format(Sys.Date(), "%d/%m/%Y"))
  }

  SoftData[SoftData$Paciente == as.integer(numPaciente), ]$Last_Update <<- format(Sys.Date(), "%d/%m/%Y")

  openxlsx::write.xlsx(paciente, patient_XLSX_file, sheetName =  sprintf("%s", numPaciente), rowNames = FALSE)
  openxlsx::write.xlsx(softwares, patient_XLSX_file, sheetName =  "Software", rowNames = FALSE, overwrite = TRUE)

  openxlsx::write.xlsx(mitoRStats, MitoRStats_path, sheetName = "MitoRStats", rowNames = FALSE)
  openxlsx::write.xlsx(SoftData, MitoRStats_path, sheetName = "SoftData", rowNames = FALSE, overwrite = TRUE)
}

