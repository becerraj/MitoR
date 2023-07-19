#' @title Update the MitoR-DataBase stats taken from HMTVAR
#' @description MitoR-DataBase takes information for each mutation from the HMTVAR DataBase. This DataBase is constantly updated.
#' We recommend to update the information of MitoR-DataBase at least once a month. With this function, the MitoR DataBase will be updated with the last
#' available information about each mutation from HMTVAR.
#' @param AddToPatients In case you want to update the information of the individual XLSX files from the patients as well, keep AddToPatients = TRUE.
#' Please be aware that this will take a while if you have many patients to update.
#' In case you want to update only one patient or just some of them, we recommend you to create the XLSX file from the VCF again. Check generateXLSX() for more information.
#' @return The DataBase XLSX file will automatically be updated. The same will happen with the patients XLSX files in case of AddToPatients = TRUE
#' @export
#' @examples
#' updateMitoRStats_HMTVAR()
#' updateMitoRStats_HMTVAR(addToPatients = TRUE)
updateMitoRStats_HMTVAR <- function(AddToPatients = FALSE){
  # Updates the MitoR DB from the HMTVAR database
  mitor_files <- sprintf("%s/mitor", Sys.getenv('R_LIBS_USER'))
  RDS_DB <- readRDS(sprintf("%s/MitoRSoftware/RDS_DB.rds", mitor_files))

  Freq_MitoR <- RDS_DB$Freq_MitoR[[1]]

  # Looks for the data given in the Freq_MitoR dataframe on the HMTVAR API
  updates <- buscarHMTVAR(Freq_MitoR)
  Freq_MitoR$clinVAR <- updates[[1]]
  Freq_MitoR$MitoMap <- updates[[2]]
  Freq_MitoR$dbSNP <- updates[[3]]
  Freq_MitoR$Omim <- updates[[4]]
  Freq_MitoR$Disease <- updates[[5]]
  Freq_MitoR$Freq_HMTVAR <- updates[[6]]

  # If it is not the first update, it deletes the last update row first
  updateRow <- data.frame(Paciente = "HMTVAR Update", Date = format(Sys.Date(), "%d/%m/%Y"), version_BWA = "-", version_GATK = "-", version_PICARD = "-", Last_Update = "-")
  if ("HMTVAR Update" %in% SoftData$Paciente) {
    SoftData <- SoftData[!SoftData$Paciente == "HMTVAR Update", ]
  }

  # Adds the update row
  SoftData <<- rbind(SoftData, updateRow)
  openxlsx::write.xlsx(mitoRStats, MitoRStats_Path, sheetName = "MitoRStats" , rowNames = FALSE)
  openxlsx::write.xlsx(SoftData, MitoRStats_Path, sheetName = "SoftData", rowNames = FALSE, overwrite = TRUE)

  if (AddToPatients == TRUE) {
    for (i in names(RDS_DB)[2:length(RDS_DB)]) {
      tryCatch(
        expr = {
          update_freq_MitoR(i)
        },
        error = function(e) {
          message(sprintf("An error occured while updating the patient number '%s'.
                          Please, try to update it individually by using the update_freq_MitoR() function.", i))

          print(e)
        },
        warning = function(w) {
          message(sprintf("An error occured while updating the patient number '%s'.
                          Please, try to update it individually by using the update_freq_MitoR() function.", i))
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

#' @title Adds or updates MitoR DB frequencies to the patient report
#' @description Takes the MitoR DB frequencies from all the previous analysis and loads them into the patient report.
#' It is important that you already have the analysis previously made. Otherwise, the function will analyze the patient from the VCF.
#' If you did not keep the XLSX report nor the VCF file, you must do the analysis again with the patient's fasta/fastq files or its BAM file.
#' @param numPaciente Patient ID.
#' @return Updated report and a message telling the path it is saved at.
#' @export
#' @examples
#' update_freq_MitoR("47286")
update_freq_MitoR <- function(numPaciente) {
  # Adds or updates the MitoR frequencies of the patient analysis
  mitor_files <- sprintf("%s/mitor", Sys.getenv('R_LIBS_USER'))
  RDS_DB <- readRDS(sprintf("%s/MitoRSoftware/RDS_DB.rds", mitor_files))

  Freq_MitoR <- RDS_DB$Freq_MitoR[[1]]

  patient_XLSX_file <- RDS_DB[[numPaciente]]["XLSX_Path"]
  # Reads the xlsx patient analysis: 1 sheet: mutations - 2 sheet: softwares
  # Add frequencies to 1. Add the update row to 2
  paciente <- openxlsx::read.xlsx(patient_XLSX_file, sheet = 1)
  softwares <- openxlsx::read.xlsx(patient_XLSX_file, sheet = 2)

  # Takes the mutations of the patient
  mutations <- RDS_DB[[numPaciente]]["Mutations"]
  toAddFreq <- strsplit(unlist(mutations, use.names = FALSE), "/")

  # Save in a new variable the MitoR DB frequencies of the mutations found in the patient
  freqMut_MitoRStats <- c()
  for (i in (1:length(toAddFreq))){
    freqMut_MitoRStats <- c(freqMut_MitoRStats, Freq_MitoR[(Freq_MitoR$POS == as.integer(toAddFreq[[i]][2])) & (Freq_MitoR$ALT == toAddFreq[[i]][3]), ]$Freq_MitoR)
  }

  if ("Freq_MitoR" %in% colnames(paciente)){ # If its an update of the frequency
    paciente$Freq_MitoR <- freqMut_MitoRStats
    softwares$Last_Upd_FreqMitoR <- format(Sys.Date(), "%d/%m/%Y")
  } else{ # If this is the first time it adds the frequency
    paciente <- cbind(paciente, Freq_MitoR = freqMut_MitoRStats)
    softwares <- cbind(softwares, Last_Upd_FreqMitoR = format(Sys.Date(), "%d/%m/%Y"))
  }

  SoftData[SoftData$Paciente == as.integer(numPaciente), ]$Last_Update <<- format(Sys.Date(), "%d/%m/%Y")

  openxlsx::write.xlsx(paciente, patient_XLSX_file, sheetName =  sprintf("%s", numPaciente), rowNames = FALSE)
  openxlsx::write.xlsx(softwares, patient_XLSX_file, sheetName =  "Software", rowNames = FALSE, overwrite = TRUE)

  openxlsx::write.xlsx(mitoRStats, MitoRStats_path, sheetName = "MitoRStats", rowNames = FALSE)
  openxlsx::write.xlsx(SoftData, MitoRStats_path, sheetName = "SoftData", rowNames = FALSE, overwrite = TRUE)

  return(sprintf("MitoR frequencies updated. Find the report at: %s", MitoRStats_path))
}
