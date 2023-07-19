#' @title Add another patient's information to Genes DB
#' @description Adds the information about the genes affected by the CNVs of a new patient
#' and keeps track of the frequency and type of variation of each gene.
#' @param CNV_calls dataframe which contains the information about the CNVs of the new patient.
#' @param id number or character that indicates which patient is.
#' @param ControlDB database of control patients.
#' @param PatientsDB database of patients to analyze or analyzed.
#' @param GenesDB database of genes affected.
#' @return Genes DB updated.
#' @import dplyr

update_DB_cnv <- function(CNV_calls, id, ControlDB, PatientsDB, GenesDB) {

  #If the patient was analyzed before, delete the old data:
  if(id %in% colnames(GenesDB)) {
    indic <- which(colnames(GenesDB) == id)
    GenesDB <- subset(GenesDB, select = -indic)
  }

  x <- c()
  e_var <- c()
  for (i in CNV_calls$name.gene){
    x <- append(x, strsplit(i, split = ","))
  }
  for (l in 1:length(x)){
    for (g in 1:length(x[[l]])){
      e_var <- append(e_var, x[[l]][g])
    }
  }

  Complete <- data.frame("Gene" = e_var)

  #Amount of duplications or delections for each gene
  all_dup <- c()
  all_del <- c()

  #Create 2 lists: one for the duplicated genes and another with the deleted
  g_dup <- unlist(CNV_calls[CNV_calls$type == "duplication", "name.gene" ])
  if (!(length(g_dup) == 0)) {
    all_dup <- c()
    q <- c()

    for (i in g_dup) {  q <- append(q, strsplit(i, split =","))  }
    for (l in 1:length(q)) {
      for (g in 1:length(q[[l]])) { all_dup <- append(all_dup, q[[l]][g]) }
    }
  }

  g_del <- unlist(CNV_calls[CNV_calls$type == "deletion", "name.gene" ])
  if (!(length(g_del) == 0)) {
    all_del <- c()
    w <- c()

    for (i in g_del) { w <- append(w, strsplit(i, split=",")) }
    for (l in 1:length(w)){
      for (g in 1:length(w[[l]])) { all_del<-append(all_del, w[[l]][g]) }
    }
  }

  #Indicate if the genes was duplicated or deleted for each patient (column):
  Complete$Paciente[Complete$Gene %in% all_dup ] <- "duplication"
  Complete$Paciente[Complete$Gene %in% all_del ] <- "deletion"


  if (max(dim(GenesDB)) == 0) { #If its the first patient in Genes DB
    GenesDB <- data.frame(Complete)
    GenesDB$Frec_dup <- 0
    GenesDB$Frec_del <- 0
    message("The first patient has been saved in GenesDB")

  } else {
    p <- as.data.frame(Complete[, c(1, 2)])
    GenesDB <- merge(GenesDB, p, by = "Gene", all = T)
    GenesDB <- select(GenesDB, -c("Frec_dup", "Frec_del", "Frec_variant") )
    message("Another patient has been saved in GenesDB")
  }

  names(GenesDB)[names(GenesDB) == "Paciente"] <- id

  GenesDB$Frec_dup <- 0
  GenesDB$Frec_del <- 0

  for (i in 1:nrow(GenesDB)) {
    for (j in 2:(ncol(GenesDB) - 2)) {
      if (!(is.na(GenesDB[i,j])) & GenesDB[i, j] == "deletion" ) {
        GenesDB$Frec_del[i] <- GenesDB$Frec_del[i] + 1
      } else if (!(is.na(GenesDB[i,j])) & GenesDB[i, j] == "duplication") {
        GenesDB$Frec_dup[i] <- GenesDB$Frec_dup[i] + 1
      }
    }
  }

  #Total amount of variations per gene:
  GenesDB$Frec_variant <- GenesDB$Frec_dup + GenesDB$Frec_del

  saveRDS(list(ControlDB, PatientsDB, GenesDB), file = paste(Sys.getenv('R_LIBS_USER'), "/mitorDB/cnvDB.RDS", sep=""))

  return(GenesDB)

}
