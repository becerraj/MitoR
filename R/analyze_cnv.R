#' @importFrom ExomeDepth getBamCounts
#' @import tibble
#' @title Detect and analyze copy number variations (CNVs) of a patient sample.
#' @description Generate an excel report with 2 sheets: in the first one there is
#' a table with information about the CNVs detected such as location, size, and hyperlinks to online data bases.
#' On the second sheet there is a visualization of the CNVs within the whole mitochondrial genome of the patient.
#' @param path_dir path of the directory of the patient that will be analyzed. It must contain R1 and R2 (fasta or fastq files).
#' @param minoverlap is a decimal number which indicates the needed percentage
#' of overlap from which it can be declare that a region has suffered a variation.
#' @param transit transition probability: Transition probability of the hidden Markov Chain from the normal copy number state to either a deletion or a duplication.
#' This parameter is taken from ExomeDepth algorithms.
#' @return a list:
#' - Patients DB (updated if the patient that has been analyzed is new).
#' - CNV calls (dataframe with information about CNVs) ordered for plotting.
#' - GenesDB (updated data base with information about the genes that has been detected as variant.
#' - CNV calls ordered for generation of the report
#' - Graph with the duplications colored in green and deletions colored in red.
#' @export
#' @examples
#' outAnalisis <- analyze_cnv(path_dir = "~/PacientesMito/46608/BM22-46608_R1.fastq.gz", transit = 0.7)
#' PatientsDB <- outAnalisis[[1]]
#' CNV_calls <- outAnalisis[[2]]
#' GenesDB <- outAnalisis[[3]]
#' graph <- outAnalisis[[5]]


analyze_cnv <- function(path_dir, minoverlap = 0.0001, transit = 1) {

  initial <- initialize_db()
  data("bedfileMito")
  ControlDB <- initial[[1]]
  PatientsDB <- initial[[2]]
  GenesDB <- initial[[3]]

  if(nrow(ControlDB) == 0){
    stop("You must have at least one patient in Control DB. Use command add_control_cnv().")
  }

  id <- basename(path_dir)

  # Avoid using the patient as a reference for itself
  if(id %in% colnames(ControlDB)) {
    indic <- which(colnames(ControlDB) == id)
    Pcounts <- ControlDB[, indic]

    #Delete the patient from ControlDB
    ControlDB <- subset(ControlDB, select = -indic)
    if (max(dim(ControlDB)) == 0 | ncol(ControlDB) == 0) { stop("You need to have at least one control patient") }

    #Add the patient to PatientsDB if its not there
    if (!(id %in% colnames(PatientsDB))) {
      if (max(dim(PatientsDB)) == 0) { #this will be the first patient in DB
        PatientsDB <- data.frame(Pcounts)
        names(PatientsDB)[names(PatientsDB) == "Pcounts"] <- id
        message("Your first patient has been saved!")

      } else { #not the first patient on the DB
        PatientsDB <- add_column(PatientsDB, Pcounts)
        names(PatientsDB)[names(PatientsDB) == "Pcounts"] <- id
        message("Another patient has been saved!")
      }
    }
  }

  if (!(id %in% colnames(PatientsDB))) { #If its a new patient---------------------------------------------------------
    message("This patients is not already on your database.The process may take a few minutes")
    bam.control <- fasta_to_bam(path_dir)
    indexBam(bam.control)
    cts <- getBamCounts(bed.frame = bed,
                        bam.files = bam.control ,
                        include.chr = F, #if set to TRUE, this function will add the string 'chr' to the chromosome names of the target BED file.
                        referenceFasta = NULL)

    Pcounts <- data.frame(Pcounts = cts[, ncol(cts)])

    #Add the new patient to the PatientsDB
    if (max(dim(PatientsDB)) == 0) { #this will be the first patient in DB
      PatientsDB <- data.frame(Pcounts)
      names(PatientsDB)[names(PatientsDB) == "Pcounts"] <- id
      message("Your first patient has been saved!")

    } else { #not the first control patient
      PatientsDB <- add_column(PatientsDB, Pcounts)
      names(PatientsDB)[names(PatientsDB) == "Pcounts"] <- id
      message("Another patient has been saved!")
    }
  }#------------------------------------------------------------------------------------------------------------------

  tryCatch(
    {
      CNV_calls <- cnv_from_counts(id, bed, ControlDB, PatientsDB, minoverlap = minoverlap, transit = transit)
      print("There are CNVs detected")
    },
    error = function(e) {
      print(paste("No CNVs were found. Try again using a higher value of transit"))
    }
  )

#}

#------------------------------------------------------------------------

  graph <- plot_cnv(CNV_calls)
  message("The graph has been succesfully generated")

  #Keep the version of CNV calls' dataframe original without reordering it
  CNV_calls_ <- CNV_calls

  #Order to generate links
  CNV_calls <- subset(CNV_calls, select = c("chromosome", "type", "size", "id", "name.gene", "nexons","BF", "reads.expected", "reads.observed", "reads.ratio"))
  colnames(CNV_calls)[1] <-"chr"

  CNV_calls <- generate_links(CNV_calls)

  wb <- generate_cnv_report(CNV_calls, graph = graph, path_dir, transit = transit, minoverlap = minoverlap)
  openxlsx::saveWorkbook(wb, paste(path_dir, "/MitoR/", id, "_CNVreport.xlsx", sep = ""),  overwrite = TRUE)
  browseURL(paste(path_dir, "/MitoR/", id, "_CNVreport.xlsx", sep = ""))
  message("Report has been succesfully generated")

  GenesDB <- update_DB_cnv(CNV_calls, id, ControlDB, PatientsDB, GenesDB)
  saveRDS(list(ControlDB, PatientsDB, GenesDB), file = paste(Sys.getenv('R_LIBS_USER'), "/mitorDB/cnvDB.RDS", sep=""))
  message("Genes DB has been succesfully updated")

  return(list(PatientsDB, CNV_calls_, GenesDB, CNV_calls, graph))
}


