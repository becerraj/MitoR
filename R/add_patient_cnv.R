#' @importFrom ExomeDepth getBamCounts
#' @import tibble
#' @title Add patient to Patients DB
#' @description Saves the counts of the lectures of each gene in a data base
#' which columns are the ids of the patients. It also changes the RDS file
#' where the databases are saved.
#' @param id number or character that indicates which patient will be added.
#' @param path_dir path of the directory of the patient that will be analyzed. It must contain R1 and R2 (fasta or fastq files).
#' @return Patients DB updated
#' @export 
#' @examples
#' PatientsDB <-  AddPatient( id="123",  path_dir ="/home/daniela/Patients/123/123.fastq")


add_patient_cnv <- function(path_dir) {
  
  initial <- initialize_db()
  data("bedfileMito")
  ControlDB <- initial[[1]]
  PatientsDB <- initial[[2]]
  GenesDB <- initial[[3]]
  
  id <- basename(path_dir)
  
  if (id %in% colnames(PatientsDB)){ #If the patients is in Patients DB
    stop("This patient is already in your database")
    
  } else if (id %in% colnames(ControlDB)) {
    indice <- which(colnames(ControlDB) == id)
    Pcounts <- ControlDB[, indice]
    message("This patient is already in your control database")
    
  } else { #If the patient is not on any DB
    message("This patient is not already in the data base. This process will take a few minutes")
    bam.control <- bam_from_fastq(path_dir)
    cts <- getBamCounts(bed.frame = bed, 
                        bam.files = bam.control ,
                        include.chr = F, #if set to TRUE, will add the string 'chr' to the chromosome names of the target BED file. 
                        referenceFasta = NULL)
    
    Pcounts <- data.frame(Pcounts=cts[, ncol(cts)])
    
  }
  
  #Save counts on PatientsDB
  if (max(dim(PatientsDB)) == 0) { #this will be the first patient in DB
    PatientsDB <- data.frame(Pcounts)
    names(PatientsDB)[names(PatientsDB) == "Pcounts"] <- id
    message("Your first patient has been saved!")
    
  } else { #not the first control patient
    PatientsDB <- add_column(PatientsDB, Pcounts)
    names(PatientsDB)[names(PatientsDB) == "Pcounts"] <- id
    message("Another patient has been saved!")
  }
  
  #Save the updated Patients DB
  #saveRDS(list(ControlDB, PatientsDB, GenesDB), file = "~/DataBases/DBs.RDS")
  saveRDS(list(ControlDB, PatientsDB, GenesDB), file = paste(Sys.getenv('R_LIBS_USER'), "/mitorDB/cnvDB.RDS", sep=""))
  
  return(PatientsDB)
}
