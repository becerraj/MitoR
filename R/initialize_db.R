#' @title Initialize Data Bases
#' @description Generation of 3 empty databases: Control(normal references),
#'  Patients (analyzed) and Genes (genes detected with variation) if it is the first
#'  time that the user uses this package. Else, it reads the databases already
#'  created.
#' @return list with 3 databases
#' @examples
#' initial <- InitializeDB()
#' ControlDB <- initial[[1]]
#' PatientsDB <- initial[[2]]
#' GenesDB <- initial[[3]]
#' @export

initialize_db <- function() {

  dir_lib <- Sys.getenv('R_LIBS_USER')

  if (!(file.exists(sprintf("%s/mitorDB", dir_lib)))) {

    dir.create(sprintf("%s/mitorDB", dir_lib))

    wb1 <- data.frame()
    wb2 <- data.frame()
    wb3 <- data.frame()

    saveRDS(list(wb1, wb2, wb3), file = paste(dir_lib, "/mitorDB/cnvDB.RDS", sep=""))
    message("Data Bases have been created. You can find them in mitorDB folder")
  }

  cnv_DBs <- readRDS(paste(Sys.getenv('R_LIBS_USER'), "/mitorDB/cnvDB.RDS", sep=""))
  ControlDB <- cnv_DBs[[1]]
  PatientsDB <- cnv_DBs[[2]]
  GenesDB <- cnv_DBs[[3]]

  return(list(ControlDB, PatientsDB, GenesDB))
}



