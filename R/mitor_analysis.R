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


mitor_analysis <- function(path_dir, cnv_minoverlap = 0.0001, cnv_transit = 1, QD_SNPS = "< 2.0", FS_SNPS = "> 60.0", MQ_SNPS = "< 40.0", MQRankSum_SNPS = "< -12.5", ReadPosRankSum_SNPS = "< -8.0", QD_INDELS = "< 2.0", FS_INDELS = "> 200.0", ReadPosRankSum_INDELS = "< -20.0") {

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

  #Obtain the sorted bam - ACA ES EL PROBLEMA CON JUAN
  path <- fasta_to_bam(path_dir)

  cts <- getBamCounts(bed.frame = bed,
                      bam.files = path ,
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
      CNV_calls <- cnv_from_counts(id, bed, ControlDB, PatientsDB, transit = transit)
      print("There are CNVs detected")
    },
    error = function(e) {
      print(paste("No CNVs were found. Try again using a higher value of transit"))
    }
  )

  graph <- plot_cnv(CNV_calls)
  CNV_calls <- subset(CNV_calls, select = c("chromosome", "type", "size", "id", "name.gene", "nexons","BF", "reads.expected", "reads.observed", "reads.ratio"))
  colnames(CNV_calls)[1] <-"chr"
  CNV_calls <- generate_links(CNV_calls)
  GenesDB <- update_DB_cnv(CNV_calls, id, ControlDB, PatientsDB, GenesDB)
  saveRDS(list(ControlDB, PatientsDB, GenesDB), file = paste(Sys.getenv('R_LIBS_USER'), "/mitorDB/cnvDB.RDS", sep=""))


  #wb_cnv <- generate_cnv_report(CNV_calls, graph, path_dir, transit, minoverlap)

  message("CNVs' analysis has been done!")

  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  #         MUTATIONS ANALYSIS
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  workbooks <- SNP_Indel_report(path_dir, QD_SNPS, FS_SNPS, MQ_SNPS, MQRankSum_SNPS, ReadPosRankSum_SNPS, QD_INDELS, FS_INDELS, ReadPosRankSum_INDELS)
  message("Mutations' analysis has been done!")

  wb1<- workbooks[[1]]
  wb2<- workbooks[[2]]
  plots_to_xlsx <- workbooks[[3]]


  # Creat complete workbook
  wb_complete <- createWorkbook()

  # Mutation's Sheet
  d1 <- read.xlsx(wb1, sheet = 1)
  addWorksheet(wb_complete, sheetName = "SNP-INDEL")
  writeData(wb_complete, sheet = "SNP-INDEL", x = d1)

  #Mutation's Plot
  addWorksheet(wb_complete, sheetName = "Plot_Mut")
  for (i in 1:length(plots_to_xlsx)) { # i is related to the X axis
    for (j in 1:length(plots_to_xlsx[[i]])) { # j s related to the Y axis
      print(plots_to_xlsx[[i]][j])
      openxlsx::insertPlot(wb_complete, 2, xy = c((i*9)-8, (j*23)-21), width = 18.5, height = 12, fileType = "png", units = "cm")
    }
  }

  #CNV's Sheet
  addWorksheet(wb_complete, sheetName = "CNV")
  writeData(wb_complete, sheet = "CNV", x = CNV_calls)

  addWorksheet(wb_complete, sheetName = "Plot_CNV")
  print(graph)
  insertPlot(wb_complete, "Plot_CNV", xy = c(1, 1), width = 25, height = 18, fileType = "png", units = "cm")

  #Software's Information Sheet
  d2 <- read.xlsx(wb2, sheet = 1)
  addWorksheet(wb_complete, sheetName = "Soft")
  writeData(wb_complete, sheet = "Soft", x = d2)

  #openxlsx::saveWorkbook(wb_complete, "wb_final.xlsx",  overwrite = TRUE)
  #browseURL("wb_final.xlsx")

  openxlsx::saveWorkbook(wb_complete, paste(dirname(path_dir),"/", id, "_REPORT.xlsx", sep = ""),  overwrite = TRUE)
  browseURL(paste(dirname(path_dir),"/MitoR/", id, "_REPORT.xlsx", sep = ""))
}






