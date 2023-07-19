#' @title Remove a Patient's Analysis from the MitoR DataBase.
#' @description After adding several patients analysis to the DataBase, you might want to remove some of them. This function allows you to remove
#' the patient from the DataBase. You also have the possibility to delete the complete analysis (except for the VCF file).
#' @param numPaciente Patient's ID or title of the patient. In case you don't remember its name, we recommend you to look for it at the second page of the XLSX DataBase.
#' In case you still don't find it, run the function using a random number as numPaciente and the system will provide you with a list of the available patients.
#' @param keepXLSX After removing the patient from the DataBase, it will also remove the individual XLSX analysis created for the patient. TRUE by default.
#' @return The patient's XLSX file containing the analysis will be removed. The patients mutations will also be removed from the DataBase if wanted.
#' @examples
#' deletePatient(37339)
#' deletePatient(37339, keepXLSX = FALSE)
deletePatient <- function(numPaciente, keepXLSX = TRUE) {
  mitor_files <- sprintf("%s/mitor", Sys.getenv('R_LIBS_USER'))
  numPaciente <- as.character(numPaciente)
  RDS_DB <- readRDS(sprintf("%s/MitoRSoftware/RDS_DB.rds", mitor_files))

  if (!RDS_DB[numPaciente] %in% names(RDS_DB)) {
    message("-------------------------------------------------------------
            The patients available in your MitoR-DataBase are: ")
    print(names(RDS_DB[2:length(RDS_DB)]))
    stop("The patient ", numPaciente, " is not in your MitoR-DataBase
         -------------------------------------------------------------")
  }

  # Removes the patient's XLSX analysis file
  if (keepXLSX == FALSE) {
    if (file.exists(as.character(RDS_DB[[numPaciente]]["XLSX_Path"]))) {
      file.remove(as.character(RDS_DB[[numPaciente]]["XLSX_Path"]))
    }}

  # Modifies the Freq_MitoR data according to the new Mutations Frequencies
  RDS_DB[1] <- deleteFreq_MitoR(numPaciente)

  # Eliminar del RDS el paciente segun el numPaciente
  # Capaz se puede hacer segun el numero de filas. Sino directamente con el numPaciente
  RDS_DB <- RDS_DB[names(RDS_DB) != (sprintf("%s", numPaciente))]
  if (length(RDS_DB) == 1) { # There's only the Freq_MitoR left
    file.remove(sprintf("%s/MitoRSoftware/RDS_DB.rds", mitor_files))
  } else {
    saveRDS(RDS_DB, "~/MitoRSoftware/RDS_DB.rds")
  }
}
