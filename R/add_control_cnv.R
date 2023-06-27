#' @importFrom ExomeDepth getBamCounts
#' @import tibble
#' @title Add patient to Control DB
#' @description Saves the counts of the lectures of each gene in a control data base 
#' which columns are the ids of the normal patients. It also changes the RDS file
#' where the databases are saved.
#' @param path_dir path of the directory of the patient that will be analyzed. It must contain R1 and R2 (fasta/q files).
#' @return Control DB updated
#' @examples ControlDB <- AddControl( id="123",  path_dir ="/home/daniela/Patients/123/123.fastq")
#' @export 
#' 
add_control_cnv <- function(path_dir) {
  initial <- initialize_db()
  data("bedfileMito")
  ControlDB <- initial[[1]]
  PatientsDB <- initial[[2]]
  GenesDB <- initial[[3]]
  
  id <- basename(path_dir)
  
  if (!(id %in% colnames(ControlDB))) { #Make sure the patient is not already in controlDB
    
    if (id %in% colnames(PatientsDB)) { #If the patients is in PatientsDB-----------------------------
      message("This patient is already in the data base. It will be copied to the control data base")
      indice <- which(colnames(PatientsDB) == id)
      Pcounts <- PatientsDB[, indice]
      
    } else { #If the patient is not on the DB------------------------------------------------------------
      message("This patient is not already in the data base. This process will take a few minutes")
      bam.control <- bam_from_fastq(path_dir)
      cts <- getBamCounts(bed.frame = bed, 
                          bam.files = bam.control ,
                          include.chr = F, # if set to TRUE, will add the string 'chr' to the chromosome names of the target BED file. 
                          referenceFasta = NULL)
      
      Pcounts <- data.frame(Pcounts = cts[, ncol(cts)])
      
    }#---------------------------------------------------------------------------
    
    #update the database
    
    if (max(dim(ControlDB)) == 0) { #this will be the first control patient
      ControlDB <- data.frame(Pcounts)
      names(ControlDB)[names(ControlDB) == "Pcounts"] <- id
      message("Your first control patient has been saved!")
      
    } else { #not the first control patient
      ControlDB <- add_column(ControlDB, Pcounts)
      names(ControlDB)[names(ControlDB) == "Pcounts"] <- id
      message("Another control patient has been saved!")
    }
    
    #saveRDS(list(ControlDB, PatientsDB, GenesDB), file = "~/DataBases/DBs.RDS")
    saveRDS(list(ControlDB, PatientsDB, GenesDB), file = paste(Sys.getenv('R_LIBS_USER'), "/mitorDB/cnvDB.RDS", sep=""))
    
  } else { stop("This patient is already on your Control Data Base") }
  
  return(ControlDB)
  
}

