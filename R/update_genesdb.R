#' @title Add another patient's information to Exons DB
#' @description Adds the information about the CNVs of exons of a new patient
#' and keeps track of the frequency and type of variation of each exon.
#' @param CNV_calls dataframe which contains the information about the CNVs of the new patient
#' @param id number or character that indicates which patient is.
#' @param e name of the RDS file  for GenesDB
#' @return Genes DB updated
#' @import dplyr

update_DB_cnv <- function(CNV_calls, id, ControlDB, PatientsDB, GenesDB) {
  x <- c()
  exones_var <- c()
  for (i in CNV_calls$name.gene){
    x <- append(x, strsplit(i, split = ","))
  }
  for (l in 1:length(x)){
    for (g in 1:length(x[[l]])){
      exones_var<-append(exones_var, x[[l]][g])
    }
  }
  
  Completo <- data.frame("Gene" = exones_var)
  
  #Amount of duplications or delections for each gene
  todos_dup <- c()
  todos_del <- c()
  
  #Create 2 lists: one for the duplicated genes and another with the deleted
  exones_duplicados <- unlist(CNV_calls[CNV_calls$type == "duplication", "name.gene" ])
  if (!(length(exones_duplicados) == 0)) {
    todos_dup <- c()
    q <- c()
    
    for (i in exones_duplicados) {  q <- append(q, strsplit(i, split =","))  }
    for (l in 1:length(q)) {
      for (g in 1:length(q[[l]])) { todos_dup <- append(todos_dup, q[[l]][g]) }
    }
  }
  
  exones_delecionados <- unlist(CNV_calls[CNV_calls$type == "deletion", "name.gene" ])
  if (!(length(exones_delecionados) == 0)) {
    todos_del <- c()
    w <- c()
    
    for (i in exones_delecionados) { w <- append(w, strsplit(i, split=",")) }
    for (l in 1:length(w)){
      for (g in 1:length(w[[l]])) { todos_del<-append(todos_del, w[[l]][g]) }
    }
  }
  
  
  #Indicate if the genes was duplicated or deleted for each patient (column): 
  Completo$Paciente[Completo$Gene %in% todos_dup ] <- "duplication"
  Completo$Paciente[Completo$Gene %in% todos_del ] <- "deletion"
  
  
  if (max(dim(GenesDB)) == 0) { #If its the first patient in Exons DB
    GenesDB <- data.frame(Completo)
    GenesDB$Frec_dup <- 0
    GenesDB$Frec_del <- 0
    message("The first patient has been saved in GenesDB")
    
  } else {
    p <- as.data.frame(Completo[, c(1, 2)])
    GenesDB <- merge(GenesDB, p, by = "Gene", all = T)
    GenesDB <- select(GenesDB, -c("Frec_dup", "Frec_del", "Frec_variant") ) ##
    message("another patient has been saved in GenesDB")
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
  
  #Total amount of variations per gene
  GenesDB$Frec_variant <- GenesDB$Frec_dup + GenesDB$Frec_del
  
  
  saveRDS(list(ControlDB, PatientsDB, GenesDB), file = "~/DataBases/DBs.RDS")
  
  return(GenesDB)
  
}