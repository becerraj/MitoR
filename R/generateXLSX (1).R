#' @title Convert VCF file to XLSX file
#' @description It generates a XLSX file with the information of the mutations available at the VCF file.
#' Creates a XLSX file from a VCF file. If you do not have one, we recommend you to see mitor::generateVCF()
#' By doing this, the genetist or the bioinformatitian is able to have the information on a more common file type.
#' Information about the mutation is also taken from HMTVAR. To use this function you must be online.
#'
#' The function has also the possibility to store the data in a database automatically generated. The database will also be a XLSX file.
#'
#' @param addToDataBase=TRUE While creating the XLSX file, the information is added to a DataBase.
#' If you want to keep add it as part of the analysis, it will be saved in the upper directory of the R1 and R2 files of the analyzed patient.
#' @return XLSX file
#' @import stringr
#' @import rjson
#' @import openxlsx
#' @export
#' @examples
#' generateXLSX("../MitoR/patient1.vcf, addToDataBase=FALSE")
#' generateXLSX("../MitoR/patient1.vcf")

#VCF_file <- "/home/juan/Github/MitochondriaAnalysis/Pacientes/37339/MitoR/37339_filtered_MitoR.vcf"
generateXLSX <- function(VCF_file) {
  numPaciente <- unlist(stringr::str_split(basename(VCF_file), "_"))[1]

  # Leemos, tomamos los parametros de filtrado y sacamos el header del VCF
  data_vcf <- analyzeVCF(VCF_file)

  # Hay que agrupar todos los parametros en un solo string asi lo agregamos en la segunda pagina del XLSX
  params <- as.character(data_vcf[2])
  patient_vcf <- unlist(data_vcf[1])

  # Preparamos el nombre de las columnas del XLSX: "CHROM" "POS" "ID" "REF" "ALT" "QUAL" "FILTER" "INFO" "FORMAT" "bar"
  col_names <- c(unlist(strsplit(patient_vcf[1], "\t")))
  col_names[1] <- "CHROM"

  # Preparamos los datos del XLSX: [1] "MT" [2] "73" [3] "." [4] "A" ...
  datos <- unlist(strsplit(patient_vcf[2:length(patient_vcf)], "\t"))

  # Unimos el nombre de las columnas y los datos en una matriz. Luego convertimos a data.frame
  matriz <- t(matrix(datos, nrow = length(col_names), ncol=(length(patient_vcf)-1)))
  filtered_DF <- as.data.frame(matriz, row.names = FALSE)
  colnames(filtered_DF) <- col_names

  # Guardamos las versiones de los Softwares usados
  versionGATK <- unlist(stringr::str_split(GATK, "[-/]"))[length(unlist(stringr::str_split(GATK, "[-/]")))-1]
  versionPICARD <- unlist(stringr::str_split(PICARD, "[/-]"))[length(unlist(stringr::str_split(PICARD, "[/-]")))-1]
  versionBWA <- "0.7.17"

  # Separamos la informacion de la columna bar: GT AD DP GQ PL
  GT <- c()
  AD <- c()
  DP <- c()
  GQ <- c()
  PL <- c()

  for (i in (1:nrow(filtered_DF))){
    dividido <- strsplit(filtered_DF$bar[i], ":")
    GT <- c(GT, dividido[[1]][1])
    AD <- c(AD, dividido[[1]][2])
    DP <- c(DP, dividido[[1]][3])
    GQ <- c(GQ, dividido[[1]][4])
    PL <- c(PL, dividido[[1]][5])
  }

  # Formateamos GT para saber si es una mutacion homocigota o heterocigota
  GT <- homHet(GT)
  Franklin <- sprintf("https://franklin.genoox.com/clinical-db/variant/snp/chrM-%s-%s-%s", filtered_DF$POS, filtered_DF$REF , filtered_DF$ALT)
  VarSome <- sprintf("https://varsome.com/variant/hg38/M%s%s%s%s%s%s?annotation-mode=germline", "%3A", filtered_DF$POS, "%3A", filtered_DF$REF , "%3A", filtered_DF$ALT)

  # Buscamos en  la API de HMTVAR datos referidos a cada mutacion
  HMTVAR <- buscarHMTVAR(filtered_DF)
  clinvar <- HMTVAR[[1]]
  dbSNP <- sprintf("https://www.ncbi.nlm.nih.gov/snp/%s", HMTVAR[[2]])
  MitoMap <- HMTVAR[[3]]
  Omim <- HMTVAR[[4]]
  Disease <- HMTVAR[[5]]
  All_freq_h <- substr(HMTVAR[[6]], start = 0, stop = 5)

  # Descartamos las columnas con info innecesaria
  drop <- names(filtered_DF) %in% c("ID", "FILTER", "FORMAT", "bar")
  filtered_DF <- filtered_DF[, !drop]

  # Unimos las columnas de la API
  filtered_DF <- cbind(filtered_DF, GT = GT, Allele_Depth = AD, Deep_Coverage = DP, GQ = GQ, PL = PL, clinVAR = clinvar, dbSNP = dbSNP, MitoMap = MitoMap, Omim = Omim, Disease = Disease, Franklin = Franklin, VarSome = VarSome, Freq_HMTVAR = All_freq_h)
  filtered_DF <- filtered_DF[c("CHROM", "POS", "REF", "ALT", "Freq_HMTVAR", "Allele_Depth",	"Deep_Coverage", "clinVAR",	"dbSNP",	"MitoMap",	"Omim",	"Disease", "Franklin", "VarSome", "INFO",	"GT",	"GQ",	"PL")]

  # Generamos la Sheet2 que contiene informacion del paciente, fecha y version de los Softwares
  filtered_DF2 <- data.frame(Patient = numPaciente, Date = format(Sys.Date(), "%d/%m/%Y"), Filter_Params = params,
                             VersionBWA = versionBWA, VersionGATK = versionGATK,
                             VersionPICARD = versionPICARD, Last_Update = format(Sys.Date(), "%d/%m/%Y"))

  # Path de las mutaciones filtradas del paciente
  XLSX_file <- paste(substr(VCF_file, start = 0, stop = (nchar(VCF_file)-4)), ".xlsx", sep = "", collapse = "")

  # Ademas de agregarse a la base de datos, agregamos al paciente la frecuencia de las mutaciones en nuestra base de datos actualizada.
  freq_MitoR <- unlist(add_to_RDS(XLSX_file)) # Genera el RDS de la base de datos. Devuelve las frecuencias de las mutaciones en MitoR

  # Agregamos la columna de frecuencias de MitoR
  filtered_DF <- cbind(filtered_DF, Freq_MitoR = freq_MitoR)
  filtered_DF <- filtered_DF[c("CHROM", "POS", "REF", "ALT", "Freq_MitoR", "Freq_HMTVAR",	"Allele_Depth",	"Deep_Coverage", "clinVAR",	"dbSNP",	"MitoMap",	"Omim",	"Disease", "Franklin", "VarSome", "INFO",	"GT",	"GQ",	"PL")]

  # Por ultimo, generamos los graficos de las mutaciones. Tiene que hacerse despues de haber creado
  # la DB porque la info proviene de ahi.
  plots_to_xlsx <- pileupPlot(numPaciente) %>% plot_xlsx_design()


  # Es necesario el parametro del file XLSX entonces tuvimos que crearlo si o si antes tambien.
  # Ahora que se modifico, hay que reescribirlo
  openxlsx::write.xlsx(filtered_DF, XLSX_file, sheetName =  sprintf("%s", numPaciente), rowNames = FALSE)
  openxlsx::write.xlsx(filtered_DF2, XLSX_file, sheetName = "Softwares", rowNames = FALSE, overwrite = TRUE)

  return(XLSX_file)
}

# Funciones dependientes de generateXLSX:
  # buscarHMTVAR()
  # homHet()
  # addToRDS()

#' @title buscarHMTVAR
#' @description Takes information from the HMTVAR DataBase about the mutations we are interested in.
#' The taken data for each mutation is: clinvar, dbSNP, MitoMap, Omim, Disease, All_freq_h.
#' @param dataFrame It must have at least a REF, POS and ALT colomn name.
#' @return lista
#' @import httr
#' @import rjson

buscarHMTVAR <- function(dataFrame){
  clinvar <- c()
  dbSNP <- c()
  MitoMap <- c()
  Omim <- c()
  Disease <- c()
  All_freq_h <- c()
  # Checking if we have internet conection to look for the HMTVAR values
  if (!httr::http_error(httr::GET("https://www.google.com"))){
    for (i in (1:nrow(dataFrame))){
      REF <- dataFrame[i, "REF"]
      POS <- dataFrame[i, "POS"]
      ALT <- dataFrame[i, "ALT"]
      datos <- system2("curl", sprintf("-X GET 'https://www.hmtvar.uniba.it/api/main/mutation/%s%s%s' -H 'accept: application/json' --insecure", REF, POS, ALT), stdout = TRUE, wait = TRUE)
      datos1 <- rjson::fromJSON(datos)

      if (is.null(datos1$CrossRef$clinvar)){
        clinvar <- c(clinvar, "-")
      }
      else{
        clinvar <- c(clinvar, sprintf("%s", datos1$CrossRef$clinvar))
      }
      if (is.null(datos1$CrossRef$dbSNP)){
        dbSNP <- c(dbSNP, "-")
      }
      else{
        dbSNP <- c(dbSNP, sprintf("%s", datos1$CrossRef$dbSNP))
      }
      if (is.null(datos1$CrossRef$mitomap_associated_disease)){
        MitoMap <- c(MitoMap, "-")
      }
      else{
        MitoMap <- c(MitoMap, sprintf("%s", datos1$CrossRef$mitomap_associated_disease))
      }
      if (is.null(datos1$CrossRef$omim)){
        Omim <- c(Omim, "-")
      }
      else{
        Omim <- c(Omim, sprintf("%s", datos1$CrossRef$omim))
      }
      if (is.null(datos1$CrossRef$pubs_disease)){
        Disease <- c(Disease, "-")
      }
      else{
        Disease <- c(Disease, sprintf("%s", datos1$CrossRef$pubs_disease))
      }
      if (is.null(datos1$Variab$all_freq_h)){
        All_freq_h <- c(All_freq_h, "-")
      }
      else{
        All_freq_h <- c(All_freq_h, sprintf("%s", datos1$Variab$all_freq_h))
      }
    }
  } else{
    clinvar <- "Offline"
    dbSNP <- "Offline"
    MitoMap <- "Offline"
    Omim <- "Offline"
    Disease <- "Offline"
    All_freq_h <- "Offline"
  }
  return(list(clinvar = clinvar, dbSNP = dbSNP, MitoMap = MitoMap, Omim = Omim, Disease = Disease, All_freq_h = All_freq_h))
}

#' @title homHet
#' @description Given the GATK parameters as 0/0, 0/1, 1/1, this function writes its meaning regarding the type of mutation.
#' @param GT List of GATK parameters as 0/0, 0/1, 1/1
#' @return The same list (if it is coming from the generateXLSX function, it is a dataFrame colomn) with the replaced meaning.
#' @examples homHet(list(0/0,0/1, 1/1))

homHet <- function(GT){
  for (i in (1:length(GT))){
    if (GT[i] == "0/0"){
      GT[i] <- "Hom Ref"
    } else if (GT[i] == "0/1"){
      GT[i] <- "Het"
    } else {
      GT[i] <- "Hom Alt"
    }
  }
  return(GT)
}

analyzeVCF<- function(VCF_file){
  patient_vcf <- as.list(readLines(VCF_file))

  filter_params <- c(patient_vcf[grep("ID=mitor_indel_filter", patient_vcf, ignore.case = TRUE)], patient_vcf[grep("ID=mitor_snp_filter", patient_vcf)])
  filter_params <- paste("Indel Params = ", stringr::str_extract(filter_params[1], '(?<=Description=\\").*(?=\\")'), "SNP Params = ", stringr::str_extract(filter_params[2], '(?<=Description=\\").*(?=\\")'), collapse = " ")

  patient_vcf <- patient_vcf[-1:-(grep("#CHROM",patient_vcf)-1)]
  return(list(patient_vcf, filter_params))
}

checkHMTVAR_online <- function(){
  if (!httr::http_error(httr::GET("https://www.google.com"))){
    mitor_files <- sprintf("%s/mitor", Sys.getenv('R_LIBS_USER'))
    if ("RDS_DB.rds" %in% list.files(mitor_files)){
      RDS_DB <- readRDS(sprintf("%s/RDS_DB.rds"), mitor_files)
      if (any(stringr::str_detect(RDS_DB[[1]]$Freq_HMTVAR , "Offline"))){
        Freq_MitoR <- RDS_DB[[1]]
        # Seleccionamos las columnas que tienen algun valor en Offline
        valuesOffline <- Freq_MitoR[Freq_MitoR$Freq_HMTVAR == "Offline", ]
        valuesToUpdate <- buscarHMTVAR(valuesOffline)
        # Agregamos las filas modificadas en el mismo lugar que estaban antes las offline
        Freq_MitoR[valuesOffline, ] <- valuesToUpdate
        # Volvemos a guardar en el RDS_DB el dataframe actualizado y lo guardamos como RDS
        RDS_DB$Freq_MitoR <- Freq_MitoR
        saveRDS(RDS_DB, sprintf("%s/MitoRSoftware/RDS_DB.rds", mitor_files))
      }
    }
  }
}

plot_xlsx_design <- function(plots_to_xlsx) {
  for (i in (1:length(plots_to_xlsx))) {

  }
}


#generateXLSX("/home/juan/Github/MitochondriaAnalysis/Pacientes/37019/MitoR/37019_filtered_MitoR.vcf")
#generateXLSX("/home/juan/Github/MitochondriaAnalysis/Pacientes/37339/MitoR/37339_filtered_MitoR.vcf")
#generateXLSX("/home/juan/Github/MitochondriaAnalysis/Pacientes/38069/MitoR/38069_filtered_MitoR.vcf")
#generateXLSX("/home/juan/Github/MitochondriaAnalysis/Pacientes/44562/MitoR/44562_filtered_MitoR.vcf")
#generateXLSX("/home/juan/Github/MitochondriaAnalysis/Pacientes/46608/MitoR/46608_filtered_MitoR.vcf")
#generateXLSX("/home/juan/Github/MitochondriaAnalysis/Pacientes/47053/MitoR/47053_filtered_MitoR.vcf")
#generateXLSX("/home/juan/Github/MitochondriaAnalysis/Pacientes/47286/MitoR/47286_filtered_MitoR.vcf")


#RDS_DB <- readRDS("/home/juan/MitoRSoftware/RDS_DB.rds")

