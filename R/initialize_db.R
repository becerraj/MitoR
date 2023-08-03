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

initialize_db <- function() {
  data("bedfileMito")

  if (!(file.exists(sprintf("%s/mitorDB/DB/cnvDB.RDS", Sys.getenv('R_LIBS_USER'))))) {

    dir.create(sprintf("%s/mitorDB/DB", Sys.getenv('R_LIBS_USER')))

    wb1 <- data.frame()
    wb2 <- data.frame("Gene" = bed[, 4])
    wb3 <- data.frame()

    saveRDS(list(wb1, wb2, wb3), file = paste(Sys.getenv('R_LIBS_USER'), "/mitorDB/DB/cnvDB.RDS", sep=""))
    message("Data Bases have been created. You can find them in mitorDB folder")
  }

  cnv_DBs <- readRDS(paste(Sys.getenv('R_LIBS_USER'), "/mitorDB/DB/cnvDB.RDS", sep=""))
  ControlDB <- cnv_DBs[[1]]
  PatientsDB <- cnv_DBs[[2]]
  GenesDB <- cnv_DBs[[3]]

  return(list(ControlDB, PatientsDB, GenesDB))
}



