#' @title Remove a patient from Patients DB or Control DB.
#' @description Removes the patient from the database indicated. Updates the RDS file where all databases are saved.
#' @param id of the patient that will be removed.
#' @param n string that indicates  the database that will be changed. It must be
#' either "Control" or "Patients".
#' @return DB updated
#' @export
#' @examples
#' PatientsDB <- RemoveAnyPatient( id= "37556", n = "Patients")
#' ControlDB <- RemoveAnyPatient( id= "37556", n = "control")


delete_patient_cnv <- function(id, n) {

  initial <- InitializeDB()
  ControlDB <- initial[[1]]
  PatientsDB <- initial[[2]]
  GenesDB <- initial[[3]]

  if(!(n == "Control" | n == "Patients")) { stop('n must be the string "Control" or "Patients" ') }

  if (n == "Control" | n == "control" | n == "CONTROL" | n == "C" ) {
    if (!(id %in% colnames(ControlDB))) { stop("This patient is not on the database indicated") }
    indice <- which(colnames(ControlDB) == id)
    ControlDB <- subset(ControlDB, select = -indice)

  } else if (n == "Patients" | n == "patients" | n == "patient" | n == "PATIENTS" | n == "P" ) {
    if (!(id %in% colnames(PatientsDB))) { stop("This patient is not on the database indicated") }
    indice <- which(colnames(PatientsDB) == id)
    PatientsDB <- subset(PatientsDB, select = -indice)
  }

  saveRDS(list(ControlDB, PatientsDB, GenesDB), file = paste(Sys.getenv('R_LIBS_USER'), "/mitorDB/cnvDB.RDS", sep=""))
  message(paste("The patient ", id, " has been removed from the database ", n, sep = ""))

  return(DB)

}
