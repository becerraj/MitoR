#' @title Analyze mutations and variations 
#' @description kjbfv
#' @param path_dir path of the directory of the patient that will be analyzed. It must contain R1 and R2 (fasta or fastq files).
#' @param minoverlap is a decimal number which indicates the needed percentage
#' of overlap from which it can be declare that a region has suffered a variation.
#' @param transit transition probability: Transition probability of the hidden Markov Chain 
#' from the normal copy number state to either a deletion or a duplication. 
#' @return a list: Patients DB (updated if the patient that has been analyzed is new),
#' CNV calls (dataframe with information about CNVs),
#' GenesDB (updated data base with information about the genes that has been detected
#' as variant)
#' @export 
#' @examples
#' outAnalisis <- MitoRAnalysis(path_dir = "~/PacientesMito/46608/BM22-46608_R1.fastq.gz")
#' PatientsDB<-outAnalisis[[1]]
#' CNV_calls<-outAnalisis[[2]]
#' all.exons<-outAnalisis[[3]]


mitor_analysis <- function(path_dir, minoverlap = 0.0001, transit = 1 ) {
  
  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  #         CNVs ANALYSIS
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  
  initial <- initialize_db()
  data("bedfileMito")
  ControlDB <- initial[[1]]
  PatientsDB <- initial[[2]]
  GenesDB <- initial[[3]]
  
  id <- basename(path_dir)
  
  # Avoid using the patient as a reference for itself
  if(id %in% colnames(ControlDB)) { 
    indice <- which(colnames(ControlDB) == id)
    Pcounts <- ControlDB[, indice]
    
    #Delete the patient from ControlDB
    ControlDB <- subset(ControlDB, select = -indice)
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
  
  if (id %in% colnames(PatientsDB)) { #If the patients is in Patients DB---------------
    message("This patients is already on your database.The process will take less than a minute")
    tryCatch(
      {
        out_counts <- cnv_from_counts(id, bed, ControlDB, PatientsDB, transit = transit)
        print("There are CNVs detected")
        CNV_calls <- out_counts[[1]]
      },
      error = function(e) {
        print(paste("No CNVs were found. Try again using a higher value of transit"))
      }
    )
  } else { #If its a new patient---------------------------------------------------------
    
    message("This patients is not already on your database.The process may take a few minutes")
    bam.control <- bam_from_fastq(path_dir)
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
      
    } else { #not the first control patient ESTO ADELANTE
      PatientsDB <- add_column(PatientsDB, Pcounts)
      names(PatientsDB)[names(PatientsDB) == "Pcounts"] <- id
      message("Another patient has been saved!")
    }
    
    tryCatch(
      {
        out_counts <- cnv_from_counts(id, bed, ControlDB, PatientsDB, transit = transit)
        print("There are CNVs detected")
        CNV_calls <- out_counts[[1]]
      },
      error = function(e) {
        print(paste("No CNVs were found. Try again using a higher value of transit"))
      }
    )
    
  }
  
  #------------------------------------------------------------------------
  
  graph <- plot_cnv(CNV_calls, bed)
  CNV_calls <- subset(CNV_calls, select = c("chromosome", "type", "size", "id", "name.gene", "nexons","BF", "reads.expected", "reads.observed", "reads.ratio"))
  colnames(CNV_calls)[1] <-"chr"
  
  CNV_calls <- generate_links(CNV_calls)
  
  GenesDB <- update_DB_cnv(CNV_calls, id, ControlDB, PatientsDB, GenesDB)
  saveRDS(list(ControlDB, PatientsDB, GenesDB), file = paste(Sys.getenv('R_LIBS_USER'), "/mitorDB/cnvDB.RDS", sep=""))
  
  message("CNVs' analysis has been done!")
  
  
  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  #         MUTATIONS ANALYSIS
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  
  #vcf_analysis <- generate_SNPs(path, keep_BAM = FALSE, QD_SNPS = "< 2.0", FS_SNPS = "> 60.0", MQ_SNPS = "< 40.0", MQRankSum_SNPS = "< -12.5", ReadPosRankSum_SNPS = "< -8.0", QD_INDELS = "< 2.0", FS_INDELS = "> 200.0", ReadPosRankSum_INDELS = "< -20.0")
  message("Mutations' analysis has been done!")
  
  
  #------------------------------------------------------------------------------------------------------
  #------------------------------------------------------------------------------------------------------
  wb <- openxlsx::createWorkbook()
  
  COLNAMES_STYLE <- createStyle(
    fontSize = 12,
    textDecoration = "bold",
    halign = "center", valign = "center", border = "TopBottom",
    borderColour = "black")
  
  
  hyperlink_style <- openxlsx::createStyle(fontColour = "#0000FF")
  
  addWorksheet(wb, paste(id, "CNVAnalisis"))
  addWorksheet(wb, paste(id, "CNVplots"))
  writeData(wb, sheet = 1, x = cnv,
            headerStyle = COLNAMES_STYLE,
            borderStyle = "dashed",
            borders = "columns", borderColour = "black")
  
  
  setColWidths(wb, sheet = 1, cols = 1:4, widths = "auto")
  setColWidths(wb, sheet = 1, cols = 6:ncol(cnv), widths = "auto")
  setColWidths(wb, sheet = 1, cols = 5, widths = 16)
  
  for(i in 1:nrow(cnv)){
    x <- cnv$Franklin[i]
    y <- cnv$Varsome[i]
    names(x) <- c("View on Franklin DB")
    names(y) <- c("View on Varsome DB")
    class(x) <- "hyperlink"
    class(y) <- "hyperlink"
    
    indiceF <- which(colnames(cnv) == "Franklin")
    indiceV <- which(colnames(cnv) == "Varsome")
    
    writeData(wb, sheet = 1, x = x, startRow= i+1, startCol = indiceF)
    writeData(wb, sheet = 1, x = y, startRow= i+1, startCol = indiceV)
    
  }
  
  #Colour cells if the size of the variation is >500
  
  for (i in 1:nrow(cnv)){
    if (cnv$size[i] > 500){
      addStyle(wb, 1, style = openxlsx::createStyle(fgFill = "thistle"), rows = i+1, cols = 1:ncol(cnv))
    }
  }
  
  #Add Plot
  print(cnv_graph)
  insertPlot(wb, 2, xy = c(1, 1), width = 25, height = 18, fileType = "png", units = "cm")
  
  
  openxlsx::saveWorkbook(wb, paste(dirname(path),"/", id, "_REPORT.xlsx", sep = ""),  overwrite = TRUE)
  browseURL(paste(dirname(path),"/MitoR/", id, "_REPORT.xlsx", sep = ""))
}
